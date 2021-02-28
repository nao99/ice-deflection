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
    private final Double value;

    /**
     * Lambda constructor
     *
     * @param value a lambda value
     * @throws IllegalArgumentException if lambda value is nullable
     */
    public Lambda(Double value) {
        if (null == value) {
            throw new IllegalArgumentException("Lambda value must be non null");
        }

        this.value = value;
    }

    public Double getValue() {
        return value;
    }
}
