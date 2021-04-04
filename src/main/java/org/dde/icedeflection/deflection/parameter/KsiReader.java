package org.dde.icedeflection.deflection.parameter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Scanner;

/**
 * KsiReader class <br>
 *
 * This class reads ksi parameters
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class KsiReader {
    /**
     * Scanner of a file with ksi parameters
     */
    private final Scanner parametersFileScanner;

    /**
     * KsiReader constructor
     *
     * @param parametersFilePath a path to file with physical parameters
     *
     * @throws IllegalArgumentException if a path is nullable
     * @throws FileNotFoundException    if a file is not valid
     */
    public KsiReader(String parametersFilePath) throws FileNotFoundException {
        if (null == parametersFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.parametersFileScanner = new Scanner(new File(parametersFilePath));
    }

    /**
     * KsiReader constructor
     *
     * @param parametersInputStream a physical parameters input stream
     * @throws IllegalArgumentException if an input stream is nullable
     */
    public KsiReader(InputStream parametersInputStream) {
        if (null == parametersInputStream) {
            throw new IllegalArgumentException("Input stream must be a non null");
        }

        this.parametersFileScanner = new Scanner(parametersInputStream);
    }

    /**
     * Reads all data from a file
     *
     * @return ksi data
     */
    public Ksi read() {
        double lowerBoundary = parametersFileScanner.nextDouble();
        double upperBoundary = parametersFileScanner.nextDouble();
        int stepsNumber = parametersFileScanner.nextInt();

        return new Ksi(lowerBoundary, upperBoundary, stepsNumber);
    }

    /**
     * Closes this reader
     */
    public void close() {
        parametersFileScanner.close();
    }
}
