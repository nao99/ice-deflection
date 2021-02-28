package org.nao99.icedeflection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.*;

/**
 * LambdaReader class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class LambdaReader {
    /**
     * Scanner of a lambda file
     */
    private final Scanner lambdaFileScanner;

    /**
     * LambdaReader constructor
     *
     * @param lambdaFilePath a path to lambda file
     *
     * @throws IllegalArgumentException if a path is nullable
     * @throws FileNotFoundException    if a file is not valid
     */
    public LambdaReader(String lambdaFilePath) throws FileNotFoundException {
        if (null == lambdaFilePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.lambdaFileScanner = new Scanner(new File(lambdaFilePath));
    }

    /**
     * LambdaReader constructor
     *
     * @param lambdasInputStream a lambdas input stream
     * @throws IllegalArgumentException if an input stream is nullable
     */
    public LambdaReader(InputStream lambdasInputStream) {
        if (null == lambdasInputStream) {
            throw new IllegalArgumentException("Input stream must be a non null");
        }

        this.lambdaFileScanner = new Scanner(lambdasInputStream);
    }

    /**
     * Reads all data from a file
     *
     * @param lambdasCount a count of lambda to read
     * @return read lambdas in object representation
     */
    public List<Lambda> read(int lambdasCount) {
        List<Lambda> lambdas = new ArrayList<>();

        int lambdasReceivedCount = 0;
        while (lambdaFileScanner.hasNext() && lambdasCount != lambdasReceivedCount) {
            double lambdaValue = lambdaFileScanner.nextDouble();
            lambdas.add(new Lambda(lambdaValue));

            lambdasReceivedCount++;
        }

        return lambdas;
    }

    /**
     * Closes this reader
     */
    public void close() {
        lambdaFileScanner.close();
    }
}
