package org.dde.icedeflection.deflection.parameter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.Scanner;

/**
 * PsiReader class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-04-04
 */
public class PsiReader {
    /**
     * Scanner of a file with psi parameters
     */
    private final Scanner parametersFileScanner;

    /**
     * PsiReader constructor
     *
     * @param parametersFilePath a path to file with physical parameters
     *
     * @throws IllegalArgumentException if a path is nullable
     * @throws FileNotFoundException    if a file is not valid
     */
    public PsiReader(String parametersFilePath) throws FileNotFoundException {
        if (null == parametersFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.parametersFileScanner = new Scanner(new File(parametersFilePath));
    }

    /**
     * PsiReader constructor
     *
     * @param parametersInputStream a physical parameters input stream
     * @throws IllegalArgumentException if an input stream is nullable
     */
    public PsiReader(InputStream parametersInputStream) {
        if (null == parametersInputStream) {
            throw new IllegalArgumentException("Input stream must be a non null");
        }

        this.parametersFileScanner = new Scanner(parametersInputStream);
    }

    /**
     * Reads all data from a file
     *
     * @return psi data
     */
    public Psi read() {
        int number = parametersFileScanner.nextInt();
        double lowerBoundary = parametersFileScanner.nextDouble();
        double upperBoundary = parametersFileScanner.nextDouble();
        int stepsNumber = parametersFileScanner.nextInt();

        return new Psi(number, lowerBoundary, upperBoundary, stepsNumber);
    }

    /**
     * Closes this reader
     */
    public void close() {
        parametersFileScanner.close();
    }
}
