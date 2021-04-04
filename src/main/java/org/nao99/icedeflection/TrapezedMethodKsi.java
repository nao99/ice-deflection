package org.nao99.icedeflection;

import java.text.DecimalFormat;
import java.util.function.Function;

/**
 * TrapezedMethod class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-22
 */
public class TrapezedMethodKsi {
    /**
     * Equation
     */
    private final Function<double[], Double> function;

    /**
     * Integration limit (a)
     */
    private final double integrationLimitLower;

    /**
     * Integration limit (b)
     */
    private final double integrationLimitUpper;

    /**
     * Steps number (n)
     */
    private final int stepsNumber;

    /**
     * TrapezedMethod constructor
     *
     * @param function              a function
     * @param integrationLimitLower a lower integration limit
     * @param integrationLimitUpper an upper integration limit
     * @param stepsNumber           a steps number
     */
    public TrapezedMethodKsi(
        Function<double[], Double> function,
        double integrationLimitLower,
        double integrationLimitUpper,
        int stepsNumber
    ) {
        this.function = function;
        this.integrationLimitLower = integrationLimitLower;
        this.integrationLimitUpper = integrationLimitUpper;
        this.stepsNumber = stepsNumber;
    }

    /**
     * Solves an {@link Function}
     *
     * @return a solution
     */
    public double solve() {
        if (integrationLimitUpper < integrationLimitLower) {
            throw new IllegalArgumentException("The upper limit of integration cannot be less than the lower");
        }

        double h = (integrationLimitUpper - integrationLimitLower) / stepsNumber;

        double x0 = integrationLimitLower;
        double xn = x0 + h;

        DecimalFormat decimalFormat = new DecimalFormat("#.######");

        double result = 0.0d;
        double index = 0.0d;

        while (xn < integrationLimitUpper) {
            result += function.apply(new double[] {xn, index});
            xn += h;

            xn = Double.parseDouble(decimalFormat.format(xn));
            index++;
        }

        result += (function.apply(new double[] {x0, 0.0d}) + function.apply(new double[] {xn, index})) / 2;
        result *= h;

        return result;
    }
}