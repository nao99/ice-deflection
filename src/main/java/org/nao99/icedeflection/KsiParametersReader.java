package org.nao99.icedeflection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Scanner;

/**
 * KsiParametersReader class <br>
 *
 * This class reads ksi parameters
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class KsiParametersReader {
    /**
     * Scanner of a file with ksi parameters
     */
    private final Scanner parametersFileScanner;

    /**
     * KsiParametersReader constructor
     *
     * @param parametersFilePath a path to file with physical parameters
     *
     * @throws IllegalArgumentException if a path is nullable
     * @throws FileNotFoundException    if a file is not valid
     */
    public KsiParametersReader(String parametersFilePath) throws FileNotFoundException {
        if (null == parametersFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.parametersFileScanner = new Scanner(new File(parametersFilePath));
    }

    /**
     * KsiParametersReader constructor
     *
     * @param parametersInputStream a physical parameters input stream
     * @throws IllegalArgumentException if an input stream is nullable
     */
    public KsiParametersReader(InputStream parametersInputStream) {
        if (null == parametersInputStream) {
            throw new IllegalArgumentException("Input stream must be a non null");
        }

        this.parametersFileScanner = new Scanner(parametersInputStream);
    }

    /**
     * Reads all data from a file
     *
     * @return read physical parameters in object representation
     */
    public KsiParameters read() {
        double lowerBoundary = parametersFileScanner.nextDouble();
        double upperBoundary = parametersFileScanner.nextDouble();
        int stepsNumber = parametersFileScanner.nextInt();

        return new KsiParameters(lowerBoundary, upperBoundary, stepsNumber);
    }

    /**
     * Closes this reader
     */
    public void close() {
        parametersFileScanner.close();
    }
}
