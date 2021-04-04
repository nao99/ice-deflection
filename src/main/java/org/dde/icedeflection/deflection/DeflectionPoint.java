package org.dde.icedeflection.deflection;

/**
 * DeflectionPoint class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-04-04
 */
public class DeflectionPoint {
    /**
     * X
     */
    private final double x;

    /**
     * Y
     */
    private final double y;

    /**
     * Z
     */
    private final double z;

    /**
     * DeflectionPoint constructor
     *
     * @param x an x coordinate
     * @param y an y coordinate
     * @param z a z coordinate
     */
    public DeflectionPoint(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return String.format("(%f, %f, %f)", x, y, z);
    }
}
