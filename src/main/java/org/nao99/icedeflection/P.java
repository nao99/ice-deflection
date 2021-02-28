package org.nao99.icedeflection;

import java.util.Map;

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
    /**
     * P function values (Map<x, y>)
     */
    private final Map<Double, Double> pMap;

    /**
     * P constructor
     *
     * @param pMap P function values
     * @throws IllegalArgumentException if p values map is nullable
     */
    public P(Map<Double, Double> pMap) {
        if (null == pMap) {
            throw new IllegalArgumentException("P values must be non null");
        }

        this.pMap = pMap;
    }

    /**
     * Gets an y value by an x value
     *
     * @param x an x value
     * @return an y value
     */
    public Double getYByX(Double x) {
        return pMap.get(x);
    }
}
