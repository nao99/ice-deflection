package org.dde.icedeflection.deflection.parameter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Scanner;

/**
 * PReader class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class PReader {
    /**
     * Scanner of a P file
     */
    private final Scanner pFileScanner;

    /**
     * PReader constructor
     *
     * @param pFilePath a path to P file
     *
     * @throws IllegalArgumentException if a path is nullable
     * @throws FileNotFoundException    if a file is not valid
     */
    public PReader(String pFilePath) throws FileNotFoundException {
        if (null == pFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.pFileScanner = new Scanner(new File(pFilePath));
    }

    /**
     * PReader constructor
     *
     * @param pValuesInputStream a P values input stream
     * @throws IllegalArgumentException if an input stream is nullable
     */
    public PReader(InputStream pValuesInputStream) {
        if (null == pValuesInputStream) {
            throw new IllegalArgumentException("Input stream must be a non null");
        }

        this.pFileScanner = new Scanner(pValuesInputStream);
    }

    /**
     * Reads all data from a file
     *
     * @return read P values in object representation
     */
    public P read() {
        double[] x = new double[1000];
        double[] y = new double[1000];

        int i = 0;
        while (pFileScanner.hasNext()) {
            String point = pFileScanner.next();
            String[] pointValues = point.split(",");

            x[i] = Double.parseDouble(pointValues[0]);
            y[i] = Double.parseDouble(pointValues[1]);

            i++;
        }

        return new P(Arrays.copyOf(x, i), Arrays.copyOf(y, i));
    }

    /**
     * Closes this reader
     */
    public void close() {
        pFileScanner.close();
    }
}
