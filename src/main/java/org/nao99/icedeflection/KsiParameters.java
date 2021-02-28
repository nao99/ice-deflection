package org.nao99.icedeflection;

import java.util.stream.Stream;

/**
 * KsiParameters class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class KsiParameters {
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
     * Ksi values
     */
    private Double[] ksiValues;

    /**
     * KsiParameters constructor
     *
     * @param lowerBoundary a lower boundary
     * @param upperBoundary an upper boundary
     * @param stepsNumber   a steps number
     */
    public KsiParameters(double lowerBoundary, double upperBoundary, int stepsNumber) {
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

    /**
     * Gets a ksi values array
     *
     * @return a ksi values array
     */
    public Double[] getKsiValues() {
        if (null != ksiValues) {
            return ksiValues;
        }

        ksiValues = Stream
            .iterate(lowerBoundary, s -> Math.round((s + step) * 100.0) / 100.0)
            .limit(stepsNumber)
            .toArray(Double[]::new);

        return ksiValues;
    }
}
