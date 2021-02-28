package org.nao99.icedeflection;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.DoubleStream;

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
    public static void main(String[] args) throws IOException {
        // 0. Entire parameters:
        int psiN = Integer.parseInt(args[0]);
        double psiY = Double.parseDouble(args[1]);

        ClassLoader classloader = Thread.currentThread().getContextClassLoader();

        InputStream physicalParametersInputStream = classloader.getResourceAsStream("physical_parameters.txt");
        InputStream p1ValuesInputStream = classloader.getResourceAsStream("p1_parameters.txt");
        InputStream p2ValuesInputStream = classloader.getResourceAsStream("p2_parameters.txt");
        InputStream lambdaEvenInputStream = classloader.getResourceAsStream("lambda_even_parameters.txt");
        InputStream lambdaOddInputStream = classloader.getResourceAsStream("lambda_odd_parameters.txt");
        InputStream ksiInputStream = classloader.getResourceAsStream("ksi_parameters.txt");

        PhysicalParametersReader parametersReader = new PhysicalParametersReader(physicalParametersInputStream);
        KsiParametersReader ksiParametersReader = new KsiParametersReader(ksiInputStream);

        PReader p1Reader = new PReader(p1ValuesInputStream);
        PReader p2Reader = new PReader(p2ValuesInputStream);

        LambdaReader lambdaEvenReader = new LambdaReader(lambdaEvenInputStream);
        LambdaReader lambdaOddReader = new LambdaReader(lambdaOddInputStream);

        PhysicalParameters physicalParameters = parametersReader.read();
        KsiParameters ksiParameters = ksiParametersReader.read();

        P P1 = p1Reader.read();
        P P2 = p2Reader.read();

        List<Lambda> lambdasEven = lambdaEvenReader.read(psiN);
        List<Lambda> lambdasOdd = lambdaOddReader.read(psiN);

        parametersReader.close();
        ksiParametersReader.close();

        p1Reader.close();
        p2Reader.close();

        lambdaEvenReader.close();
        lambdaOddReader.close();

        int yStepsNumber = (int) ((Y_BOUNDARY_RIGHT - Y_BOUNDARY_LEFT) / psiY);

        // 3. Calculate A and B values for odd and even cases (they are needed to calculate C and M matrices)
        double[] BArrayEven = lambdasEven.stream().flatMapToDouble(l -> DoubleStream.of(cos(l.getValue()) / cosh(l.getValue()))).toArray();
        double[] AArrayEven = Arrays.stream(BArrayEven).map(B -> sqrt(1 / (1 + B * B))).toArray();

        double[] BArrayOdd = lambdasOdd.stream().flatMapToDouble(l -> DoubleStream.of(sin(l.getValue()) / sinh(l.getValue()))).toArray();
        double[] AArrayOdd = Arrays.stream(BArrayOdd).map(B -> sqrt(1 / (1 - B * B))).toArray();

        // 4. Calculate psi (vogues) (ψc[j], ψs[j], where c - even, s - odd)
        //RealMatrix psiMatrixEven = new Array2DRowRealMatrix(psiN, yStepsNumber + 1);
        //RealMatrix psiMatrixOdd = new Array2DRowRealMatrix(psiN, yStepsNumber + 1);

        List<Map<Double, Double>> psiMatrixEven = new ArrayList<>();
        List<Map<Double, Double>> psiMatrixOdd = new ArrayList<>();

        DecimalFormat decimalFormat = new DecimalFormat("####.######");

        for (int i = 0; i < psiN; i++) {
            double lambdaEven = lambdasEven.get(i).getValue();
            double lambdaOdd = lambdasOdd.get(i).getValue();

            double AEven = AArrayEven[i];
            double AOdd = AArrayOdd[i];

            double BEven = BArrayEven[i];
            double BOdd = BArrayOdd[i];

            Map<Double, Double> psiEvenMap = new HashMap<>();
            Map<Double, Double> psiOddMap = new HashMap<>();

            for (int j = 0; j < yStepsNumber + 1; j++) {
                double stepValue = Y_BOUNDARY_LEFT + psiY * j;
                stepValue = Double.parseDouble(decimalFormat.format(stepValue));

                double psiEven = AEven * (cos(lambdaEven * stepValue) - BEven * cosh(lambdaEven * stepValue));
                psiEvenMap.put(stepValue, psiEven);

                double psiOdd = AOdd * (sin(lambdaOdd * stepValue) - BOdd * sinh(lambdaOdd * stepValue));
                psiOddMap.put(stepValue, psiOdd);
            }

            psiMatrixEven.add(psiEvenMap);
            psiMatrixOdd.add(psiOddMap);
        }

        // 5. Calculate C matrices for even and odd cases (Cc[j], Cs[j], where c - even, s - odd)
        RealMatrix CMatrixEven = new Array2DRowRealMatrix(psiN, psiN);
        RealMatrix CMatrixOdd = new Array2DRowRealMatrix(psiN, psiN);

        for (int n = 0; n < psiN; n++) {
            double lambdaEvenN = lambdasEven.get(n).getValue();
            double lambdaOddN = lambdasOdd.get(n).getValue();

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
                double lambdaEvenM = lambdasEven.get(m).getValue();
                double lambdaOddM = lambdasOdd.get(m).getValue();

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

        RealMatrix UJMatrixEven = new Array2DRowRealMatrix(j, ksiParameters.getStepsNumber());
        RealMatrix UJMatrixOdd = new Array2DRowRealMatrix(j, ksiParameters.getStepsNumber());

        RealMatrix FIJMatrixEven = new Array2DRowRealMatrix(j, psiN);
        RealMatrix FIJMatrixOdd = new Array2DRowRealMatrix(j, psiN);

        for (int jIndex = 0; jIndex < j; jIndex++) {
            for (int psiIndex = 0; psiIndex < psiN; psiIndex++) {
                double lambdaEven = lambdasEven.get(psiIndex).getValue();
                double lambdaOdd = lambdasOdd.get(psiIndex).getValue();

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
            for (double ksi : ksiParameters.getKsiValues()) {
                double uEven = sqrt(ksi * ksi + pow(PI * jIndex, 2));
                double uOdd = sqrt(ksi * ksi + pow(PI * jIndex - PI / 2, 2));

                UJMatrixEven.addToEntry(jIndex, ksiIndex, uEven);
                UJMatrixOdd.addToEntry(jIndex, ksiIndex, uOdd);

                ksiIndex++;
            }
        }

        int ksiIndex = 0;
        for (double ksi : ksiParameters.getKsiValues()) {
            RealMatrix MMatrixEven = new Array2DRowRealMatrix(psiN, psiN);
            RealMatrix MMatrixOdd = new Array2DRowRealMatrix(psiN, psiN);

            double ksiSquared = ksi * ksi;

            for (int m = 0; m < psiN; m++) {
                for (int n = 0; n < psiN; n++) {
                    double MEven = FIJMatrixEven.getEntry(0, m) * FIJMatrixEven.getEntry(0, n)
                        / (2 * ksi * tanh(ksi * physicalParameters.getHDimensionless()));

                    if (0.0 == ksi) {
                        MEven = 0.0;
                    }

                    double MOdd = 0.0;

                    for (int jIndex = 1; jIndex < j; jIndex++) {
                        double ujEven = UJMatrixEven.getEntry(jIndex, ksiIndex);
                        double ujOdd = UJMatrixOdd.getEntry(jIndex, ksiIndex);

                        MEven += FIJMatrixEven.getEntry(jIndex, m) * FIJMatrixEven.getEntry(jIndex, n) / ujEven * tanh(ujEven * physicalParameters.getHDimensionless());
                        MOdd += FIJMatrixOdd.getEntry(jIndex, m) * FIJMatrixOdd.getEntry(jIndex, n) / ujOdd * tanh(ujOdd * physicalParameters.getHDimensionless());
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

        for (double ksi : ksiParameters.getKsiValues()) {
            RealMatrix DMatrixEven = new Array2DRowRealMatrix(psiN, psiN);
            RealMatrix DMatrixOdd = new Array2DRowRealMatrix(psiN, psiN);

            double ksiFourth = pow(ksi, 4);

            for (int i = 0; i < psiN; i++) {
                double DEven = pow(lambdasEven.get(i).getValue(), 4) + ksiFourth;
                double DOdd = pow(lambdasOdd.get(i).getValue(), 4) + ksiFourth;

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
        for (double ksi : ksiParameters.getKsiValues()) {
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
        for (double ksi : ksiParameters.getKsiValues()) {
            RealMatrix QMatrixEven = QMatricesEven.get(m);
            RealMatrix QMatrixOdd = QMatricesOdd.get(m);

            RealMatrix MMatrixEven = MMatricesEven.get(m);
            RealMatrix MMatrixOdd = MMatricesOdd.get(m);

            RealMatrix IMatrix = new Array2DRowRealMatrix(psiN, psiN);
            for (int n = 0; n < psiN; n++) {
                IMatrix.addToEntry(n, n, 1 - physicalParameters.getAlpha() * physicalParameters.getHDimensionless() * physicalParameters.getFr() * physicalParameters.getFr() * ksi * ksi);
            }

            RealMatrix QMatrixEvenMultiplied = QMatrixEven.scalarMultiply(physicalParameters.getBeta());
            RealMatrix QMatrixOddMultiplied = QMatrixOdd.scalarMultiply(physicalParameters.getBeta());

            RealMatrix MMatrixEvenMultiplied = MMatrixEven.scalarMultiply(physicalParameters.getHDimensionless() * physicalParameters.getFr() * physicalParameters.getFr());
            RealMatrix MMatrixOddMultiplied = MMatrixOdd.scalarMultiply(physicalParameters.getHDimensionless() * physicalParameters.getFr() * physicalParameters.getFr());

            RealMatrix RMatrixEven = IMatrix.add(QMatrixEvenMultiplied).subtract(MMatrixEvenMultiplied);
            RealMatrix RMatrixOdd = IMatrix.add(QMatrixOddMultiplied).subtract(MMatrixOddMultiplied);

            RMatricesEven.add(RMatrixEven);
            RMatricesOdd.add(RMatrixOdd);

            m++;
        }

        // 10. Calculate P1F arrays
        double[] P1FArray = new double[ksiParameters.getStepsNumber()];

        final double integrationLimitLower = 0.0;
        final double integrationLimitUpper = 0.2;
        final int stepsNumber = 320; // will be moved to params

        for (int ksiI = 0; ksiI < ksiParameters.getStepsNumber(); ksiI++) {
            double ksi = ksiParameters.getKsiValues()[ksiI];
            TrapezedMethod trapezedMethod = new TrapezedMethod(
                (x) -> P1.getYByX(x) * cos(ksi * x),
                integrationLimitLower,
                integrationLimitUpper,
                stepsNumber
            );

            P1FArray[ksiI] = sqrt(2 / PI) * trapezedMethod.solve();
        }

        // 11. Calculate Pm* array
        double[] PMEvenArray = new double[psiN];
        double[] PMOddArray = new double[psiN];

        RealMatrix PMMatrixEven = new Array2DRowRealMatrix(psiN, 1);
        RealMatrix PMMatrixOdd = new Array2DRowRealMatrix(psiN, 1);

        for (int psiI = 0; psiI < psiN; psiI++) {
            int psiIFinal = psiI;

            TrapezedMethod trapezedMethodEven = new TrapezedMethod(
                (x) -> P2.getYByX(x) * psiMatrixEven.get(psiIFinal).get(x),
                -1.0,
                1.0,
                160
            );

            TrapezedMethod trapezedMethodOdd = new TrapezedMethod(
                (x) -> P2.getYByX(x) * psiMatrixOdd.get(psiIFinal).get(x),
                -1.0,
                1.0,
                160
            );

            PMMatrixEven.addToEntry(psiI, 0, trapezedMethodEven.solve());
            PMMatrixOdd.addToEntry(psiI, 0, trapezedMethodOdd.solve());

            PMEvenArray[psiIFinal] = trapezedMethodEven.solve();
            PMOddArray[psiIFinal] = trapezedMethodOdd.solve();
        }

        // 12. Calculate P vectors for even and odd cases (Pm(ksi))
        List<RealMatrix> PEvenArray = new ArrayList<>();
        List<RealMatrix> POddArray = new ArrayList<>();

        for (int ksiI = 0; ksiI < ksiParameters.getStepsNumber(); ksiI++) {
            double P1F = P1FArray[ksiI];

            RealMatrix PmKsiMatrixEven = PMMatrixEven.scalarMultiply(P1F);
            RealMatrix PmKsiMatrixOdd = PMMatrixOdd.scalarMultiply(P1F);

            PEvenArray.add(PmKsiMatrixEven);
            POddArray.add(PmKsiMatrixOdd);
        }

        // 13. Calculate al, ar for even and odd cases
        List<RealMatrix> aRMatricesEven = new ArrayList<>();
        List<RealMatrix> aRMatricesOdd = new ArrayList<>();
        List<RealMatrix> aIMatricesEven = new ArrayList<>();
        List<RealMatrix> aIMatricesOdd = new ArrayList<>();

        int n = 0;
        for (double ksi : ksiParameters.getKsiValues()) {
            // 12.1 Calculate RFinal matrix for even and odd cases
            RealMatrix RMatrixEven = RMatricesEven.get(n);
            RealMatrix RMatrixOdd = RMatricesOdd.get(n);

            RealMatrix PMatrixEven = PEvenArray.get(n);
            RealMatrix PMatrixOdd = POddArray.get(n);

            RealMatrix QMatrixEven = QMatricesEven.get(n);
            RealMatrix QMatrixOdd = QMatricesOdd.get(n);

            RealMatrix RMatrixEvenInverse = MatrixUtils.inverse(RMatrixEven);
            RealMatrix RMatrixOddInverse = MatrixUtils.inverse(RMatrixOdd);

            RealMatrix QMatrixEvenSquared = QMatrixEven.multiply(QMatrixEven);
            RealMatrix QMatrixOddSquared = QMatrixOdd.multiply(QMatrixOdd);

            RealMatrix QSquaredRInverseMatrixEven = QMatrixEvenSquared.multiply(RMatrixEvenInverse);
            RealMatrix QSquaredRInverseMatrixOdd = QMatrixOddSquared.multiply(RMatrixOddInverse);

            RealMatrix RFinalMatrixEven = RMatrixEven.add(QSquaredRInverseMatrixEven.scalarMultiply(physicalParameters.getBeta() * physicalParameters.getBeta() * ksi * ksi * physicalParameters.getEpsilon() * physicalParameters.getE()));
            RealMatrix RFinalMatrixOdd = RMatrixOdd.add(QSquaredRInverseMatrixOdd.scalarMultiply(physicalParameters.getBeta() * physicalParameters.getBeta() * ksi * ksi * physicalParameters.getEpsilon() * physicalParameters.getE()));

            // 12.2 Calculate RFinal inverse negative matrix for even and odd cases
            RealMatrix RFinalMatrixInverseNegativeEven = MatrixUtils.inverse(RFinalMatrixEven).scalarMultiply(- 1.0);
            RealMatrix RFinalMatrixInverseNegativeOdd = MatrixUtils.inverse(RFinalMatrixOdd).scalarMultiply(- 1.0);

            // 12.3 Calculate ar for even and odd cases
            RealMatrix aRMatrixEven = RFinalMatrixInverseNegativeEven.multiply(PMatrixEven);
            RealMatrix aRMatrixOdd = RFinalMatrixInverseNegativeOdd.multiply(PMatrixOdd);

            // 12.4 Calculate al for even and odd cases
            RealMatrix aIMatrixEven = RMatrixEvenInverse.multiply(QMatrixEven).multiply(aRMatrixEven).scalarMultiply(physicalParameters.getBeta() * ksi * physicalParameters.getEpsilon());
            RealMatrix aIMatrixOdd = RMatrixOddInverse.multiply(QMatrixOdd).multiply(aRMatrixOdd).scalarMultiply(physicalParameters.getBeta() * ksi * physicalParameters.getEpsilon());

            aRMatricesEven.add(aRMatrixEven);
            aRMatricesOdd.add(aRMatrixOdd);
            aIMatricesEven.add(aIMatrixEven);
            aIMatricesOdd.add(aIMatrixOdd);

            n++;
        }

//        // 13. Calculate W matrices for even and odd cases
//        BufferedWriter wWriter = new BufferedWriter(new FileWriter("/home/glen/Projects/ice_deflection/data/w_values.txt"));
//
//        for (Double x : P1.keySet()) {
//            for (Double y : P2.keySet()) {
//                double sum = 0;
//                for (int psiI = 0; psiI < psiN; psiI++) {
//                    double psiEven = psiMatrixEven.get(psiI).get(y);
//                    double psiOdd = psiMatrixOdd.get(psiI).get(y);
//
//                    int finalPsiI = psiI;
//                    TrapezedMethod wEven = new TrapezedMethod(
//                        (xn) -> aRMatricesEven.get(0).getEntry(finalPsiI,0) * cos(0 * xn) - aIMatricesEven.get(0).getEntry(finalPsiI,0) * sin(0 * xn),
//                        0,
//                        ksiStepsNumber,
//                        ksiStepsNumber
//                    );
//
//                    TrapezedMethod wOdd = new TrapezedMethod(
//                        (xn) -> aRMatricesOdd.get(0).getEntry(finalPsiI,0) * cos(0 * xn) - aIMatricesOdd.get(0).getEntry(finalPsiI,0) * sin(0 * xn),
//                        0,
//                        ksiStepsNumber,
//                        ksiStepsNumber
//                    );
//
//                    sum += psiEven * wEven.solve() + psiOdd * wOdd.solve();
//                }
//
//                sum *= sqrt(2 / PI);
//
//                String xyzValue = String.format("%4.6f,%4.6f,%4.6f\n", x, y, sum);
//                wWriter.write(xyzValue);
//            }
//        }
//
//        wWriter.close();
    }
}
