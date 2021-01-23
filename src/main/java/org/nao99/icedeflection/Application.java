package org.nao99.icedeflection;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

import static java.lang.Math.*;

/**
 * Application class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-01-23
 */
public class Application {
    private static final int Y_BOUNDARY_LEFT  = -1;
    private static final int Y_BOUNDARY_RIGHT = 1;

    /**
     * The application entry point
     *
     * @param args arguments
     */
    public static void main(String[] args) {
        // 1. Entire the main parameters:
        //  - voguesN = vogues count (voguesN > 0),
        //  - x = step 1 (-1 < x < 1)
        //  - y = step 2 (-1 < y < 1)
        //  ...

        int voguesN = Integer.parseInt(args[0]);
        if (0 >= voguesN || 6 < voguesN) { // 6 boundary exists only because of we have 6 calculated lambdas
            throw new IllegalArgumentException(String.format("Vogues count must be positive, but %d given", voguesN));
        }

        // -1 < x < 1
        double x = Double.parseDouble(args[1]);
        if (-1 >= x || x >= 1) {
            throw new IllegalArgumentException(String.format("Incorrect x: %f (-1 < x < 1 excepted)", x));
        }

        // -1 < y < 1
        double y = Double.parseDouble(args[2]);
        if (Y_BOUNDARY_LEFT >= y || y >= Y_BOUNDARY_RIGHT) {
            throw new IllegalArgumentException(String.format("Incorrect y: %f (-1 < y < 1 excepted)", y));
        }

        y = 0.0125;
        int stepsNumber = (int) ((Y_BOUNDARY_RIGHT - Y_BOUNDARY_LEFT) / y);

        // 2. Calculate odd and even lambdas (λc[j], λs[j], where c - even, s - odd)
        // Since we have already calculated them, we won't do it again. Just copy it
        // We don't need recalculate them because they are always constants

        double[] lambdasEven = {
            2.365020372431352,
            5.497803919000836,
            8.63937982869974,
            11.78097245102023,
            14.92256510455163,
            18.06415775814131,
        };

        double[] lambdasOdd = {
            3.927378719118806,
            7.068584195523235,
            10.21017612552063,
            13.35176877775915,
            16.49336143134642,
            19.63495408493621,
        };

        // 3. Calculate A and B values for odd and even cases (it need to calculate C and M matrices in the future)
        // Bc[j] = cos(λc[j]) / cosh(λc[j]), Ac[j] = sqrt(1 / (1 + Bc[j]^2))
        // Bs[j] = sin(λs[j]) / sinh(λs[j]), As[j] = sqrt(1 / (1 - Bs[j]^2))

        double[] BValuesEven = Arrays.stream(lambdasEven).map(l -> cos(l) / cosh(l)).toArray();
        double[] AValuesEven = Arrays.stream(BValuesEven).map(B -> sqrt(1 / (1 + B * B))).toArray();

        double[] BValuesOdd = Arrays.stream(lambdasOdd).map(l -> sin(l) / sinh(l)).toArray();
        double[] AValuesOdd = Arrays.stream(BValuesOdd).map(B -> sqrt(1 / (1 - B * B))).toArray();

        // 4. Calculate vogues (ψc[j], ψs[j], where c - even, s - odd)

        RealMatrix voguesEven = new Array2DRowRealMatrix(voguesN, stepsNumber);
        RealMatrix voguesOdd = new Array2DRowRealMatrix(voguesN, stepsNumber);

        for (int j = 0; j < voguesN; j++) {
            double lambdaEven = lambdasEven[j];
            double lambdaOdd = lambdasOdd[j];

            double AEven = AValuesEven[j];
            double AOdd = AValuesOdd[j];

            double BEven = BValuesEven[j];
            double BOdd = BValuesOdd[j];

            for (int stepNumber = 0; stepNumber < stepsNumber; stepNumber++) {
                double stepValue = Y_BOUNDARY_LEFT + y * stepNumber;

                double vogueEven = AEven * (cos(lambdaEven * stepValue) - BEven * cosh(lambdaEven * stepValue));
                voguesEven.addToEntry(j, stepNumber, vogueEven);

                double vogueOdd = AOdd * (sin(lambdaOdd * stepValue) - BOdd * sinh(lambdaOdd * stepValue));
                voguesOdd.addToEntry(j, stepNumber, vogueOdd);
            }
        }

        // 5. Calculate C matrices for even and odd cases (Cc[j], Cs[j], where c - even, s - odd)

        RealMatrix CValuesEven = new Array2DRowRealMatrix(voguesN, voguesN);
        RealMatrix CValuesOdd = new Array2DRowRealMatrix(voguesN, voguesN);

        for (int n = 0; n < voguesN; n++) {
            double lambdaEvenN = lambdasEven[n];
            double lambdaOddN = lambdasOdd[n];

            double AEvenN = AValuesEven[n];
            double AOddN = AValuesOdd[n];

            double lambdaEvenNSquare = lambdaEvenN * lambdaEvenN;
            double lambdaOddNSquare = lambdaOddN * lambdaOddN;

            double AEvenNSquare = AEvenN * AEvenN;
            double AOddNSquare = AOddN * AOddN;

            double lambdaEvenNCos = cos(lambdaEvenN);
            double lambdaOddNCos = cos(lambdaOddN);

            double lambdaEvenNSin = sin(lambdaEvenN);
            double lambdaOddNSin = sin(lambdaOddN);

            for (int m = 0; m < voguesN; m++) {
                double lambdaEvenM = lambdasEven[m];
                double lambdaOddM = lambdasOdd[m];

                double lambdaEvenMSquare = lambdaEvenM * lambdaEvenM;
                double lambdaOddMSquare = lambdaOddM * lambdaOddM;

                double AEvenM = AValuesEven[m];
                double AOddM = AValuesOdd[m];

                double lambdaEvenMCos = cos(lambdaEvenM);
                double lambdaOddMCos = cos(lambdaOddM);

                double lambdaEvenMSin = sin(lambdaEvenM);
                double lambdaOddMSin = sin(lambdaOddM);

                double CEven;
                double COdd;

                // TODO: We have a little error (find an fix it)
                if (m == n) {
                    double lambdaEvenNCosSquare = lambdaEvenNCos * lambdaEvenNCos;
                    double lambdaEvenNCoshSquare = cosh(lambdaEvenN) * cosh(lambdaEvenN);

                    double lambdaOddNSinSquare = lambdaOddNSin * lambdaOddNSin;
                    double lambdaOddNSinhSquare = sinh(lambdaOddN) * sinh(lambdaOddN);

                    CEven = lambdaEvenNSquare * AEvenNSquare;
                    CEven *= lambdaEvenNCosSquare / lambdaEvenNCoshSquare - 1
                        - (2 * lambdaEvenNSin * lambdaEvenNCos) / lambdaEvenN;

                    COdd = lambdaOddNSquare * AOddNSquare;
                    COdd *= (- lambdaOddNSinSquare / lambdaOddNSinhSquare) - 1
                        - (2 * lambdaOddNSin * lambdaOddNCos) / lambdaOddN;
                } else {
                    CEven = 8 * lambdaEvenNSquare * lambdaEvenMSquare * AEvenN * AEvenM;
                    CEven /= lambdaEvenNSquare * lambdaEvenNSquare - lambdaEvenMSquare * lambdaEvenMSquare;
                    CEven *= lambdaEvenM * lambdaEvenNCos * lambdaEvenMSin - lambdaEvenN * lambdaOddMCos * lambdaEvenNSin;

                    COdd = 8 * lambdaOddNSquare * lambdaOddMSquare * AOddN * AOddM;
                    COdd /= lambdaOddNSquare * lambdaOddNSquare - lambdaOddMSquare * lambdaOddMSquare;
                    COdd *= lambdaOddN * lambdaOddNCos * lambdaOddMSin - lambdaOddM * lambdaEvenMCos * lambdaOddNSin;
                }

                CValuesEven.addToEntry(n, m, CEven);
                CValuesOdd.addToEntry(n, m, COdd);
            }
        }
    }
}
