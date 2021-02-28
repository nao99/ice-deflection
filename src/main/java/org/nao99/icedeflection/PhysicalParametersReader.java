package org.nao99.icedeflection;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * PhysicalParametersReader class <br>
 *
 * This class reads initial physical parameters
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class PhysicalParametersReader {
    /**
     * Scanner of a file with physical parameters
     */
    private final Scanner parametersFileScanner;

    /**
     * PhysicalParametersReader constructor
     *
     * @param parametersFilePath a path to file with physical parameters
     *
     * @throws IllegalArgumentException if a path is nullable
     * @throws FileNotFoundException    if a file is not valid
     */
    public PhysicalParametersReader(String parametersFilePath) throws FileNotFoundException {
        if (null == parametersFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.parametersFileScanner = new Scanner(new File(parametersFilePath));
    }

    /**
     * Reads all data from a file
     *
     * @return read physical parameters in object representation
     */
    public PhysicalParameters read() {
        double hi = parametersFileScanner.nextDouble();
        double L = parametersFileScanner.nextDouble();
        double H = parametersFileScanner.nextDouble();
        double E = parametersFileScanner.nextDouble();
        double nu = parametersFileScanner.nextDouble();
        double tau = parametersFileScanner.nextDouble();
        double pi = parametersFileScanner.nextDouble();
        double pl = parametersFileScanner.nextDouble();
        double U = parametersFileScanner.nextDouble();

        return new PhysicalParameters(hi, L, H, E, nu, tau, pi, pl, U);
    }

    /**
     * Closes this reader
     */
    public void close() {
        parametersFileScanner.close();
    }
}
