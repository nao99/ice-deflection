package org.dde.icedeflection.deflection.parameter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
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
public class PhysicReader {
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
    public PhysicReader(String parametersFilePath) throws FileNotFoundException {
        if (null == parametersFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.parametersFileScanner = new Scanner(new File(parametersFilePath));
    }

    /**
     * PhysicalParametersReader constructor
     *
     * @param parametersInputStream a physical parameters input stream
     * @throws IllegalArgumentException if an input stream is nullable
     */
    public PhysicReader(InputStream parametersInputStream) {
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
    public Physic read() {
        double hi = parametersFileScanner.nextDouble();
        double L = parametersFileScanner.nextDouble();
        double H = parametersFileScanner.nextDouble();
        double E = parametersFileScanner.nextDouble();
        double nu = parametersFileScanner.nextDouble();
        double tau = parametersFileScanner.nextDouble();
        double pi = parametersFileScanner.nextDouble();
        double pl = parametersFileScanner.nextDouble();
        double U = parametersFileScanner.nextDouble();

        return new Physic(hi, L, H, E, nu, tau, pi, pl, U);
    }

    /**
     * Closes this reader
     */
    public void close() {
        parametersFileScanner.close();
    }
}
