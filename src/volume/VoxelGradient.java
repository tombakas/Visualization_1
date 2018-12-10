/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class VoxelGradient {

    public float x, y, z;
    public float mag;
    
    public VoxelGradient() {
        x = y = z = mag = 0.0f;
    }
    
    public VoxelGradient(float gx, float gy, float gz) {
        x = gx;
        y = gy;
        z = gz;
        mag = (float) Math.sqrt(x*x + y*y + z*z);
    }

    public VoxelGradient mult(float t) {
        return new VoxelGradient(x * t, y * t, z * t);
    }

    public VoxelGradient add(VoxelGradient add) {
        return new VoxelGradient(x + add.x, y + add.y, z + add.z);
    }
}
