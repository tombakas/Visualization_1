/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    static String renderFunction = "slicer";

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public static void setRenderFunction(String newRenderFunction) {
        renderFunction = newRenderFunction;
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for stering the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }


    double getTrilinearVoxel(double[] coord) {
        if (coord[0] < 0 || Math.ceil(coord[0]) >= volume.getDimX() ||
            coord[1] < 0 || Math.ceil(coord[1]) >= volume.getDimY() ||
            coord[2] < 0 || Math.ceil(coord[2]) >= volume.getDimZ())
        {
            return getVoxel(coord);
        }

        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        double alpha, alpha_expression;
        double beta, beta_expression;
        double gamma, gamma_expression;

        alpha_expression = (x - Math.floor(x))/(Math.ceil(x) - Math.floor(x));
        beta_expression = (y - Math.floor(y))/(Math.ceil(y) - Math.floor(y));
        gamma_expression = (z - Math.floor(z))/(Math.ceil(z) - Math.floor(z));

        alpha = (!Double.isNaN((alpha_expression))) ? alpha_expression : 0.0;
        beta = (!Double.isNaN(beta_expression)) ? beta_expression : 0.0;
        gamma = (!Double.isNaN(gamma_expression)) ? gamma_expression : 0.0;

        double x0 = volume.getVoxel((int) Math.floor(x), (int) Math.floor(y), (int) Math.floor(z));
        double x1 = volume.getVoxel((int) Math.ceil(x), (int) Math.floor(y), (int) Math.floor(z));
        double x2 = volume.getVoxel((int) Math.floor(x), (int) Math.floor(y), (int) Math.ceil(z));
        double x3 = volume.getVoxel((int) Math.ceil(x), (int) Math.floor(y), (int) Math.ceil(z));
        double x4 = volume.getVoxel((int) Math.floor(x), (int) Math.ceil(y), (int) Math.floor(z));
        double x5 = volume.getVoxel((int) Math.ceil(x), (int) Math.ceil(y), (int) Math.floor(z));
        double x6 = volume.getVoxel((int) Math.floor(x), (int) Math.ceil(y), (int) Math.ceil(z));
        double x7 = volume.getVoxel((int) Math.ceil(x), (int) Math.ceil(y), (int) Math.ceil(z));

        double voxelVal = (1 - alpha) * (1 - beta) * (1 - gamma) * x0 + alpha * (1 - beta) * (1 - gamma) * x1 +
                (1 - alpha) * beta * (1 - gamma) * x2 + alpha * beta * (1 - gamma) * x3 +
                (1 - alpha) * (1 - beta) * gamma * x4 + alpha * (1 - beta) * gamma * x5 +
                (1 - alpha) * beta * gamma * x6 + alpha * beta * gamma*x7;
        return voxelVal;
    }


    public int [] getEntryExit(int i, int j, int maxDist, double [] uVec, double [] vVec, double [] viewVec, double[] volumeCenter) {
        int [] entryExit = new int[2];
        double val;

        // image is square
        int imageCenter = image.getWidth() / 2;
        short threshold = 0;
        short step = 1;

        double[] pixelCoord = new double[3];

        for (int k = 0; k < maxDist; k += step) {
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

            val = getVoxel(pixelCoord);
            if (val > threshold) {
                if (k > step) {
                    entryExit[0] = k - step;
                } else {
                    entryExit[0] = k;
                }
                break;
            }
        }

        for (int k = maxDist; k > 0; k -= step) {
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                    + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                    + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                    + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

            val = getVoxel(pixelCoord);
            if (val > threshold) {
                entryExit[1] = k + step + 1;
                if (k != maxDist) {
                    entryExit[0] = k - step;
                } else {
                    entryExit[0] = k;
                }
                break;
            }
        }
        return entryExit;
    }

    void mip(double[] viewMatrix) {
        mip(viewMatrix, false);
    }

    TFColor computeCompositeColor(TFColor [] rayColors) {
        int k = rayColors.length;

        double red = 0;
        double green = 0;
        double blue = 0;
        double alpha = 1;

        double alphaProduct;

        for (int m=0; m<k; m++) {
//            alphaProduct = 1;

//            for (int n=m+1; n<m; n++) {
//                alphaProduct *= (1 - rayColors[n].a);
//            }

            red   =  (1 - rayColors[m].a) * red   + rayColors[m].r * rayColors[m].a;
            green =  (1 - rayColors[m].a) * green + rayColors[m].g * rayColors[m].a;
            blue  =  (1 - rayColors[m].a) * blue  + rayColors[m].b * rayColors[m].a;
            alpha =  (1 - rayColors[m].a) * alpha + rayColors[m].a;
        }

        return new TFColor(red, green, blue,  alpha);
    }

    void mip(double[] viewMatrix, boolean useCompositing) {

        // clear image
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        double maxVox;
        int maxDist = (int)Math.sqrt(Math.pow(volume.getDimX(), 2) + Math.pow(volume.getDimZ(), 2) + Math.pow(volume.getDimY(), 2));

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                maxVox = 0;
                double val = 0;
                int [] entryExit;

                int step = 1;
                TFColor [] rayColors = new TFColor[(int)Math.ceil(256 / step)];

                // entryExit = getEntryExit(i, j, maxDist, uVec, vVec, viewVec, volumeCenter);

                for (int k = 0; k < 256; k+=step) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (k - volumeCenter[0]) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (k - volumeCenter[1]) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (k - volumeCenter[2]) + volumeCenter[2];

                    val = getTrilinearVoxel(pixelCoord);

                    if (useCompositing) {
                        rayColors[k / step] = tFunc.getColor((int)val);
                    } else {
                        if (val > maxVox) {
                            maxVox = val;
                        } else {
                            val = maxVox;
                        }

                    }
                }

                if (!useCompositing) {
                    voxelColor.r = val / max;
                    voxelColor.g = voxelColor.r;
                    voxelColor.b = voxelColor.r;
                    voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                } else {
                    voxelColor = computeCompositeColor(rayColors);
                }

                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                // BufferedImage expects a pixel color packed as ARGB in an int
            }
        }
    }


    void slicer(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();


        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                double val = getTrilinearVoxel(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color


                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }


    @Override
    public void visualize(GL2 gl) {

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        switch (renderFunction) {
            case "slicer": slicer(viewMatrix);
                break;
            case "mip": mip(viewMatrix);
                break;
            case "mip_comp": mip(viewMatrix, true);
                break;
        }

        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
