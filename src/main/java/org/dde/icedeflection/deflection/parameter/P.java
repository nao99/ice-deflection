package org.dde.icedeflection.deflection.parameter;

import java.util.Arrays;

import static java.lang.Math.abs;

/**
 * P class <br>
 *
 * This class contains P function as a map <br>
 * P function - is an outside load to an ice place
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class P {
    private static final double PRECISION = 0.000001;

    /**
     * X array
     */
    private final double[] x;

    /**
     * Y array
     */
    private final double[] y;

    /**
     * P constructor
     *
     * @param x an x array
     * @param y an y array
     *
     * @throws IllegalArgumentException x.length != y.length or x.length == 0
     */
    public P(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("X and Y dimensions must be equal");
        }

        if (x.length == 0) {
            throw new IllegalArgumentException("X array must not be empty");
        }

        this.x = x;
        this.y = y;
    }

    public double[] getX() {
        return x;
    }

    public double[] getY() {
        return y;
    }

    /**
     * Gets min X
     *
     * @return a min x
     */
    public double getXMin() {
        return Arrays.stream(x)
            .min()
            .getAsDouble();
    }

    /**
     * Gets min X
     *
     * @return a min x
     */
    public double getXMax() {
        return Arrays.stream(x)
            .max()
            .getAsDouble();
    }

    /**
     * Gets size of P
     *
     * @return a P size
     */
    public int size() {
        return x.length;
    }

    /**
     * Gets an y value by x
     *
     * @param xValue an x value
     *
     * @return an y value
     * @throws IllegalArgumentException if x not found in array
     */
    public double getYValue(double xValue) {
        int xIndex = -1;
        for (int i = 0; i < x.length; i++) {
            if (abs(x[i] - xValue) <= PRECISION) {
                xIndex = i;
                break;
            }
        }

        if (xIndex == -1) {
            throw new IllegalArgumentException(String.format("Unable to find \"%f\" x", xValue));
        }

        return y[xIndex];
    }
}
