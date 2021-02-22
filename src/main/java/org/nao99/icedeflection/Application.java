package org.nao99.icedeflection;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import static java.lang.Math.*;

/**
 * Application class <br>
 *
 * All variables called as common letters specially for
 * more accurate match with the task itself
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-01-23
 */
public class Application {
    private static final double G = 9.8;
    private static final int Y_BOUNDARY_LEFT  = -1;
    private static final int Y_BOUNDARY_RIGHT = 1;

    /**
     * The application entry point
     *
     * @param args arguments
     */
    public static void main(String[] args) throws FileNotFoundException {
        // 0. Entire physical parameters:
        double hi = 0.1;
        double L = 10.0;
        double H = 2.0;
        double E = 4.2 * pow(10, 9);
        double tau = 0.1;
        double nu = 0.3;
        double D = (E * pow(hi, 3)) / (12 * (1 - pow(nu, 2)));
        double pi = 917.0;
        double p = 1000;
        double U = 7.0;
        double p0 = 1000;

        double M = pi * hi;
        double h = H / L;
        double beta = D / (p * G * pow(L, 4));
        double epsilon = (tau * U) / L;
        double alpha = M / (p * L);
        double Fr = U / sqrt(G * H);

        // 0.1. Read P1 values (-0.2 < P1(x) < 0.2 | 640 steps)
        Scanner P1Scanner = new Scanner(new File("/home/glen/Projects/ice_deflection/data/p1_values.txt"));
        double[] P1 = new double[641];
        for (int i = 0; i < P1.length; i++) {
            P1[i] = P1Scanner.nextDouble();
        }

        // 0.2. Read P2 values (-1 < P2(y) < 1 | 160 steps)
        Scanner P2Scanner = new Scanner(new File("/home/glen/Projects/ice_deflection/data/p2_values.txt"));
        double[] P2 = new double[161];
        for (int i = 0; i < P2.length; i++) {
            P2[i] = P2Scanner.nextDouble();
        }

        // 1. Entire the main parameters:
        int psiN = Integer.parseInt(args[0]);
        double psiY = Double.parseDouble(args[1]);

        int yStepsNumber = (int) ((Y_BOUNDARY_RIGHT - Y_BOUNDARY_LEFT) / psiY);

        int ksiStepsNumber = 1000;
        double ksiStep = 100.0 / ksiStepsNumber;

        double[] ksiArray = new double[ksiStepsNumber];
        for (int i = 0; i < ksiStepsNumber; i++) {
            ksiArray[i] = i * ksiStep;
        }

        // 2. Calculate odd and even lambdas (λc[j], λs[j], where c - even, s - odd)
        double[] lambdaArrayEven = {
            2.365020372431352,
            5.497803919000836,
            8.63937982869974,
            11.78097245102023,
            14.92256510455163,
            18.06415775814131,
        };

        double[] lambdaArrayOdd = {
            3.927378719118806,
            7.068584195523235,
            10.21017612552063,
            13.35176877775915,
            16.49336143134642,
            19.63495408493621,
        };

        // 3. Calculate A and B values for odd and even cases (they are needed to calculate C and M matrices)
        double[] BArrayEven = Arrays.stream(lambdaArrayEven).map(l -> cos(l) / cosh(l)).toArray();
        double[] AArrayEven = Arrays.stream(BArrayEven).map(B -> sqrt(1 / (1 + B * B))).toArray();

        double[] BArrayOdd = Arrays.stream(lambdaArrayOdd).map(l -> sin(l) / sinh(l)).toArray();
        double[] AArrayOdd = Arrays.stream(BArrayOdd).map(B -> sqrt(1 / (1 - B * B))).toArray();

        // 4. Calculate psi (vogues) (ψc[j], ψs[j], where c - even, s - odd)
        RealMatrix psiMatrixEven = new Array2DRowRealMatrix(psiN, yStepsNumber);
        RealMatrix psiMatrixOdd = new Array2DRowRealMatrix(psiN, yStepsNumber);

        for (int i = 0; i < psiN; i++) {
            double lambdaEven = lambdaArrayEven[i];
            double lambdaOdd = lambdaArrayOdd[i];

            double AEven = AArrayEven[i];
            double AOdd = AArrayOdd[i];

            double BEven = BArrayEven[i];
            double BOdd = BArrayOdd[i];

            for (int j = 0; j < yStepsNumber; j++) {
                double stepValue = Y_BOUNDARY_LEFT + psiY * j;

                double psiEven = AEven * (cos(lambdaEven * stepValue) - BEven * cosh(lambdaEven * stepValue));
                psiMatrixEven.addToEntry(i, j, psiEven);

                double psiOdd = AOdd * (sin(lambdaOdd * stepValue) - BOdd * sinh(lambdaOdd * stepValue));
                psiMatrixOdd.addToEntry(i, j, psiOdd);
            }
        }

        // 5. Calculate C matrices for even and odd cases (Cc[j], Cs[j], where c - even, s - odd)
        RealMatrix CMatrixEven = new Array2DRowRealMatrix(psiN, psiN);
        RealMatrix CMatrixOdd = new Array2DRowRealMatrix(psiN, psiN);

        for (int n = 0; n < psiN; n++) {
            double lambdaEvenN = lambdaArrayEven[n];
            double lambdaOddN = lambdaArrayOdd[n];

            double AEvenN = AArrayEven[n];
            double AOddN = AArrayOdd[n];

            double lambdaEvenNSquare = lambdaEvenN * lambdaEvenN;
            double lambdaOddNSquare = lambdaOddN * lambdaOddN;

            double AEvenNSquare = AEvenN * AEvenN;
            double AOddNSquare = AOddN * AOddN;

            double lambdaEvenNCos = cos(lambdaEvenN);
            double lambdaOddNCos = cos(lambdaOddN);

            double lambdaEvenNSin = sin(lambdaEvenN);
            double lambdaOddNSin = sin(lambdaOddN);

            for (int m = 0; m < psiN; m++) {
                double lambdaEvenM = lambdaArrayEven[m];
                double lambdaOddM = lambdaArrayOdd[m];

                double lambdaEvenMSquare = lambdaEvenM * lambdaEvenM;
                double lambdaOddMSquare = lambdaOddM * lambdaOddM;

                double AEvenM = AArrayEven[m];
                double AOddM = AArrayOdd[m];

                double lambdaEvenMCos = cos(lambdaEvenM);
                double lambdaOddMCos = cos(lambdaOddM);

                double lambdaEvenMSin = sin(lambdaEvenM);
                double lambdaOddMSin = sin(lambdaOddM);

                double CEven;
                double COdd;

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

                CMatrixEven.addToEntry(n, m, CEven);
                CMatrixOdd.addToEntry(n, m, COdd);
            }
        }

        // 6. Calculate M matrices for even and odd cases (Mc[j], Ms[j], where c - even, s - odd)
        List<RealMatrix> MMatricesEven = new ArrayList<>();
        List<RealMatrix> MMatricesOdd = new ArrayList<>();
        int j = 10; // j - manually chosen parameter

        RealMatrix UJMatrixEven = new Array2DRowRealMatrix(j, ksiStepsNumber);
        RealMatrix UJMatrixOdd = new Array2DRowRealMatrix(j, ksiStepsNumber);

        RealMatrix FIJMatrixEven = new Array2DRowRealMatrix(j, psiN);
        RealMatrix FIJMatrixOdd = new Array2DRowRealMatrix(j, psiN);

        for (int jIndex = 0; jIndex < j; jIndex++) {
            for (int psiIndex = 0; psiIndex < psiN; psiIndex++) {
                double lambdaEven = lambdaArrayEven[psiIndex];
                double lambdaOdd = lambdaArrayOdd[psiIndex];

                double AEven = AArrayEven[psiIndex];
                double AOdd = AArrayOdd[psiIndex];

                double FIjnEven = pow(-1, jIndex) * 4 * pow(lambdaEven, 3) * AEven * sin(lambdaEven)
                    / (pow(lambdaEven, 4) - pow(PI * jIndex, 4));

                double FIjnOdd = pow(-1, jIndex) * -4 * pow(lambdaOdd, 3) * AOdd * cos(lambdaOdd)
                    / (pow(lambdaOdd, 4) - pow(PI * jIndex - PI / 2, 4));

                FIJMatrixEven.addToEntry(jIndex, psiIndex, FIjnEven);
                FIJMatrixOdd.addToEntry(jIndex, psiIndex, FIjnOdd);
            }

            int ksiIndex = 0;
            for (double ksi : ksiArray) {
                double uEven = sqrt(ksi * ksi + pow(PI * jIndex, 2));
                double uOdd = sqrt(ksi * ksi + pow(PI * jIndex - PI / 2, 2));

                UJMatrixEven.addToEntry(jIndex, ksiIndex, uEven);
                UJMatrixOdd.addToEntry(jIndex, ksiIndex, uOdd);

                ksiIndex++;
            }
        }

        int ksiIndex = 0;
        for (double ksi : ksiArray) {
            RealMatrix MMatrixEven = new Array2DRowRealMatrix(psiN, psiN);
            RealMatrix MMatrixOdd = new Array2DRowRealMatrix(psiN, psiN);

            double ksiSquared = ksi * ksi;

            for (int m = 0; m < psiN; m++) {
                for (int n = 0; n < psiN; n++) {
                    double MEven = FIJMatrixEven.getEntry(0, m) * FIJMatrixEven.getEntry(0, n)
                        / (2 * ksi * tanh(ksi * h));

                    if (0.0 == ksi) {
                        MEven = 0.0;
                    }

                    double MOdd = 0.0;

                    for (int jIndex = 1; jIndex < j; jIndex++) {
                        double ujEven = UJMatrixEven.getEntry(jIndex, ksiIndex);
                        double ujOdd = UJMatrixOdd.getEntry(jIndex, ksiIndex);

                        MEven += FIJMatrixEven.getEntry(jIndex, m) * FIJMatrixEven.getEntry(jIndex, n) / ujEven * tanh(ujEven * h);
                        MOdd += FIJMatrixOdd.getEntry(jIndex, m) * FIJMatrixOdd.getEntry(jIndex, n) / ujOdd * tanh(ujOdd * h);
                    }

                    MEven *= ksiSquared;
                    MOdd *= ksiSquared;

                    MMatrixEven.addToEntry(m, n, MEven);
                    MMatrixOdd.addToEntry(m, n, MOdd);
                }
            }

            MMatricesEven.add(MMatrixEven);
            MMatricesOdd.add(MMatrixOdd);

            ksiIndex++;
        }

        // 7. Calculate D matrices for even and odd cases (Dc[j], Ds[j], where c - even, s - odd)
        List<RealMatrix> DMatricesEven = new ArrayList<>();
        List<RealMatrix> DMatricesOdd = new ArrayList<>();

        for (double ksi : ksiArray) {
            RealMatrix DMatrixEven = new Array2DRowRealMatrix(psiN, psiN);
            RealMatrix DMatrixOdd = new Array2DRowRealMatrix(psiN, psiN);

            double ksiFourth = pow(ksi, 4);

            for (int i = 0; i < psiN; i++) {
                double DEven = pow(lambdaArrayEven[i], 4) + ksiFourth;
                double DOdd = pow(lambdaArrayOdd[i], 4) + ksiFourth;

                DMatrixEven.addToEntry(i, i, DEven);
                DMatrixOdd.addToEntry(i, i, DOdd);
            }

            DMatricesEven.add(DMatrixEven);
            DMatricesOdd.add(DMatrixOdd);
        }

        // 8. Calculate Q matrices for even and odd cases (Qc[j], Qs[j], where c - even, s - odd)
        List<RealMatrix> QMatricesEven = new ArrayList<>();
        List<RealMatrix> QMatricesOdd = new ArrayList<>();

        int i = 0; // will be removed after code refractored and optimized
        for (double ksi : ksiArray) {
            double ksiSquared = ksi * ksi;

            RealMatrix DMatrixEven = DMatricesEven.get(i);
            RealMatrix DMatrixOdd = DMatricesOdd.get(i);

            RealMatrix KsiCMatrixEven = CMatrixEven.scalarMultiply(2 * ksiSquared);
            RealMatrix KsiCMatrixOdd = CMatrixOdd.scalarMultiply(2 * ksiSquared);

            RealMatrix QMatrixEven = DMatrixEven.subtract(KsiCMatrixEven);
            RealMatrix QMatrixOdd = DMatrixOdd.subtract(KsiCMatrixOdd);

            QMatricesEven.add(QMatrixEven);
            QMatricesOdd.add(QMatrixOdd);

            i++;
        }

        // 9. Calculate R matrices for even and odd cases
        List<RealMatrix> RMatricesEven = new ArrayList<>();
        List<RealMatrix> RMatricesOdd = new ArrayList<>();

        int m = 0;
        for (double ksi : ksiArray) {
            RealMatrix QMatrixEven = QMatricesEven.get(m);
            RealMatrix QMatrixOdd = QMatricesOdd.get(m);

            RealMatrix MMatrixEven = MMatricesEven.get(m);
            RealMatrix MMatrixOdd = MMatricesOdd.get(m);

            RealMatrix IMatrix = new Array2DRowRealMatrix(psiN, psiN);
            for (int n = 0; n < psiN; n++) {
                IMatrix.addToEntry(n, n, 1 - alpha * h * Fr * Fr * ksi * ksi);
            }

            RealMatrix QMatrixEvenMultiplied = QMatrixEven.scalarMultiply(beta);
            RealMatrix QMatrixOddMultiplied = QMatrixOdd.scalarMultiply(beta);

            RealMatrix MMatrixEvenMultiplied = MMatrixEven.scalarMultiply(h * Fr * Fr);
            RealMatrix MMatrixOddMultiplied = MMatrixOdd.scalarMultiply(h * Fr * Fr);

            RealMatrix RMatrixEven = IMatrix.add(QMatrixEvenMultiplied).subtract(MMatrixEvenMultiplied);
            RealMatrix RMatrixOdd = IMatrix.add(QMatrixOddMultiplied).subtract(MMatrixOddMultiplied);

            RMatricesEven.add(RMatrixEven);
            RMatricesOdd.add(RMatrixOdd);

            m++;
        }

        // 10. Calculate P1F arrays
        // TODO: CHECK THAT THIS METHOD CALCULATES AS WE NEED
        double[] P1FArray = new double[ksiStepsNumber];

        final double integrationLimitLower = 0.0;
        final double integrationLimitUpper = 0.2;
        final int stepsNumber = 320; // will be moved to params

        for (int ksiI = 0; ksiI < ksiStepsNumber; ksiI++) {
            double ksi = ksiArray[ksiI];
            TrapezedMethod trapezedMethod = new TrapezedMethod(
                (x) -> P1[stepsNumber] * cos(ksi * x),
                integrationLimitLower,
                integrationLimitUpper,
                stepsNumber
            );

            P1FArray[ksiI] = sqrt(2 / PI) * trapezedMethod.solve();
        }

        // 11. Calculate Pm* array
        double[] PMArray = new double[psiN];

        // 12. Calculate al, ar for even and odd cases
        int n = 0;
        for (double ksi : ksiArray) {
            // 12.1 Calculate RFinal matrix for even and odd cases
            RealMatrix RMatrixEven = RMatricesEven.get(n);
            RealMatrix RMatrixOdd = RMatricesOdd.get(n);

            RealMatrix QMatrixEven = QMatricesEven.get(n);
            RealMatrix QMatrixOdd = QMatricesOdd.get(n);

            RealMatrix RMatrixEvenInverse = MatrixUtils.inverse(RMatrixEven);
            RealMatrix RMatrixOddInverse = MatrixUtils.inverse(RMatrixOdd);

            RealMatrix QMatrixEvenSquared = QMatrixEven.multiply(QMatrixEven);
            RealMatrix QMatrixOddSquared = QMatrixOdd.multiply(QMatrixOdd);

            RealMatrix QSquaredRInverseMatrixEven = QMatrixEvenSquared.multiply(RMatrixEvenInverse);
            RealMatrix QSquaredRInverseMatrixOdd = QMatrixOddSquared.multiply(RMatrixOddInverse);

            RealMatrix RFinalMatrixEven = RMatrixEven.add(QSquaredRInverseMatrixEven.scalarMultiply(beta * beta * ksi * ksi * epsilon * epsilon));
            RealMatrix RFinalMatrixOdd = RMatrixOdd.add(QSquaredRInverseMatrixOdd.scalarMultiply(beta * beta * ksi * ksi * epsilon * epsilon));

            // 12.2 Calculate RFinal inverse negative matrix for even and odd cases
            RealMatrix RFinalMatrixInverseNegativeEven = MatrixUtils.inverse(RFinalMatrixEven).scalarMultiply(- 1.0);
            RealMatrix RFinalMatrixInverseNegativeOdd = MatrixUtils.inverse(RFinalMatrixOdd).scalarMultiply(- 1.0);

            // 12.3 Calculate ar for even and odd cases
            // 12.4 Calculate al for even and odd cases

            n++;
        }

        // 11. Calculate W matrices for even and odd cases

        String a = "";
    }
}
