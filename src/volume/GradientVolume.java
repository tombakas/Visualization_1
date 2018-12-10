/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
//        System.out.printf("%d %d %d\n", x, y, z);
        return data[x + dimX * (y + dimY * z)];
    }

    public double getTrilinearGradient(double [] coord) {
        if (coord[0] < 0 || Math.ceil(coord[0]) >= dimX ||
            coord[1] < 0 || Math.ceil(coord[1]) >= dimY ||
            coord[2] < 0 || Math.ceil(coord[2]) >= dimZ)
        {
            return 0;
        }

        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        double x_d = (x - Math.floor(x)) / (Math.ceil(x) - Math.floor(x));
        double y_d = (y - Math.floor(y)) / (Math.ceil(y) - Math.floor(y));
        double z_d = (z - Math.floor(z)) / (Math.ceil(z) - Math.floor(z));

        if (Double.isNaN(x_d)) { x_d = 0; }
        if (Double.isNaN(y_d)) { y_d = 0; }
        if (Double.isNaN(z_d)) { z_d = 0; }

        VoxelGradient c_000 = getGradient((int)Math.floor(x), (int)Math.floor(y), (int)Math.floor(z));
        VoxelGradient c_001 = getGradient((int)Math.floor(x), (int)Math.ceil(y),  (int)Math.floor(z));
        VoxelGradient c_010 = getGradient((int)Math.floor(x), (int)Math.floor(y), (int)Math.ceil(z));
        VoxelGradient c_011 = getGradient((int)Math.floor(x), (int)Math.ceil(y),  (int)Math.ceil(z));
        VoxelGradient c_100 = getGradient((int)Math.ceil(x),  (int)Math.floor(y), (int)Math.floor(z));
        VoxelGradient c_101 = getGradient((int)Math.ceil(x),  (int)Math.ceil(y),  (int)Math.floor(z));
        VoxelGradient c_110 = getGradient((int)Math.ceil(x),  (int)Math.floor(y), (int)Math.ceil(z));
        VoxelGradient c_111 = getGradient((int)Math.ceil(x),  (int)Math.ceil(y),  (int)Math.ceil(z));

        VoxelGradient c_00 = c_000.mult(1 - (float)x_d).add(c_100.mult((float)x_d));
        VoxelGradient c_01 = c_001.mult(1 - (float)x_d).add(c_101.mult((float)x_d));
        VoxelGradient c_10 = c_010.mult(1 - (float)x_d).add(c_110.mult((float)x_d));
        VoxelGradient c_11 = c_011.mult(1 - (float)x_d).add(c_111.mult((float)x_d));

        VoxelGradient c_0 = c_00.mult(1 - (float)y_d).add(c_10.mult((float)y_d));
        VoxelGradient c_1 = c_01.mult(1 - (float)y_d).add(c_11.mult((float)y_d));

        VoxelGradient c = c_0.mult(1 - (float)z_d).add(c_1.mult((float)z_d));

        return c.mag;
    }

    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    private void compute() {

        // this just initializes all gradients to the vector (0,0,0)
        for (int i=0; i<dimX; i++) {
            for (int j=0; j<dimY; j++) {
                for (int k=0; k<dimZ; k++) {
                    float x, y, z;

                    x = 0.5f * (volume.getVoxel(i+1, j, k) - volume.getVoxel(i-1, j, k));
                    y = 0.5f * (volume.getVoxel(i, j+1, k) - volume.getVoxel(i, j-1, k));
                    z = 0.5f * (volume.getVoxel(i, j, k+1) - volume.getVoxel(i, j, k-1));

                    setGradient(i, j , k, new VoxelGradient(x, y ,z));
                }
            }
        }
                
    }
    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
