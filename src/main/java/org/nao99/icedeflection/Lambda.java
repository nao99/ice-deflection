package org.nao99.icedeflection;

/**
 * Lambda class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class Lambda {
    /**
     * Lambda value
     */
    private final double value;

    /**
     * Lambda constructor
     *
     * @param value a lambda value
     */
    public Lambda(double value) {
        this.value = value;
    }

    public double getValue() {
        return value;
    }
}
