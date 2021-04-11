package org.dde.icedeflection.deflection;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

/**
 * DeflectionsWriter class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-04-11
 */
public class DeflectionsWriter {
    private static final String FORMAT = "%4.6f,%4.6f,%4.6f\n";

    /**
     * Writer
     */
    private final Writer writer;

    /**
     * DeflectionsWriter constructor
     *
     * @param filePath a file path
     * @throws IOException if an I/O error occurs
     */
    public DeflectionsWriter(String filePath) throws IOException {
        if (null == filePath) {
            throw new IllegalArgumentException("File path must be a non null");
        }

        this.writer = new BufferedWriter(new FileWriter(filePath));
    }

    /**
     * Writes {@link DeflectionPoint}s to a file
     *
     * @param points a deflection points list
     * @throws IOException if an I/O error occurs
     */
    public void write(List<DeflectionPoint> points) throws IOException {
        for (DeflectionPoint point : points) {
            String pointStr = String.format(FORMAT, point.getX(), point.getY(), point.getZ());
            writer.write(pointStr);
        }
    }

    /**
     * Closes writer
     *
     * @throws IOException if an I/O error occurs
     */
    public void close() throws IOException {
        writer.close();
    }
}
