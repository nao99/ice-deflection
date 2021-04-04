package org.dde.icedeflection.deflection.parameter;

/**
 * Ksi class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class Ksi {
    /**
     * Lower boundary
     */
    private final double lowerBoundary;

    /**
     * Upper boundary
     */
    private final double upperBoundary;

    /**
     * Steps number
     */
    private final int stepsNumber;

    /**
     * Step
     */
    private final double step;

    /**
     * Ksi constructor
     *
     * @param lowerBoundary a lower boundary
     * @param upperBoundary an upper boundary
     * @param stepsNumber   a steps number
     */
    public Ksi(double lowerBoundary, double upperBoundary, int stepsNumber) {
        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;
        this.stepsNumber = stepsNumber;
        this.step = (upperBoundary - lowerBoundary) / stepsNumber;
    }

    public double getLowerBoundary() {
        return lowerBoundary;
    }

    public double getUpperBoundary() {
        return upperBoundary;
    }

    public int getStepsNumber() {
        return stepsNumber;
    }

    public double getStep() {
        return step;
    }
}
