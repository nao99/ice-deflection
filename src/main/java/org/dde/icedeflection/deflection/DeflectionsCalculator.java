package org.dde.icedeflection.deflection;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.dde.icedeflection.deflection.parameter.Ksi;
import org.dde.icedeflection.deflection.parameter.P;
import org.dde.icedeflection.deflection.parameter.Physic;
import org.dde.icedeflection.deflection.parameter.Psi;
import org.dde.icedeflection.integral.TrapezedMethod;

import java.text.DecimalFormat;
import java.util.*;

import static java.lang.Math.*;

/**
 * DeflectionsCalculator class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-04-04
 */
public class DeflectionsCalculator {
    /**
     * Calculates deflections
     *
     * @param lambdasEven a lambdas even array
     * @param lambdasOdd  a lambdas odd array
     * @param psi         a psi data
     * @param ksi         a ksi data
     * @param physic      a physic data
     * @param P1          a p1 data
     * @param P2          a p2 data
     *
     * @return deflection points
     */
    public List<DeflectionPoint> calculate(
        double[] lambdasEven,
        double[] lambdasOdd,
        Psi psi,
        Ksi ksi,
        Physic physic,
        P P1,
        P P2,
        double[] x
    ) {
        double[] BEven = calculateBEven(lambdasEven);
        double[] BOdd = calculateBOdd(lambdasOdd);

        double[] AEven = calculateAEven(BEven);
        double[] AOdd = calculateAOdd(BOdd);

        double[] ksiValues = new double[ksi.getStepsNumber()];
        for (int i = 0; i < ksi.getStepsNumber(); i++) {
            double value = round((ksi.getLowerBoundary() + ksi.getStep() * i) * 100.0d) / 100.0d;
            ksiValues[i] = value;
        }

        // 2. Calculate psi (vogues) (ψc[j], ψs[j], where c - even, s - odd)
        List<Map<Double, Double>> psiMatrixEven = new ArrayList<>();
        List<Map<Double, Double>> psiMatrixOdd = new ArrayList<>();

        DecimalFormat decimalFormat = new DecimalFormat("####.######");
        for (int i = 0; i < psi.getNumber(); i++) {
            double lambdaEven = lambdasEven[i];
            double lambdaOdd = lambdasOdd[i];

            double AEvenValue = AEven[i];
            double AOddValue = AOdd[i];

            double BEvenValue = BEven[i];
            double BOddValue = BOdd[i];

            Map<Double, Double> psiEvenMap = new HashMap<>();
            Map<Double, Double> psiOddMap = new HashMap<>();

            for (int j = 0; j < psi.getStepsNumber() + 1; j++) {
                double stepValue = psi.getLowerBoundary() + psi.getStep() * j;
                stepValue = Double.parseDouble(decimalFormat.format(stepValue));

                double psiEven = AEvenValue * (cos(lambdaEven * stepValue) - BEvenValue * cosh(lambdaEven * stepValue));
                psiEvenMap.put(stepValue, psiEven);

                double psiOdd = AOddValue * (sin(lambdaOdd * stepValue) - BOddValue * sinh(lambdaOdd * stepValue));
                psiOddMap.put(stepValue, psiOdd);
            }

            psiMatrixEven.add(psiEvenMap);
            psiMatrixOdd.add(psiOddMap);
        }

        // 3. Calculate C matrices for even and odd cases (Cc[j], Cs[j], where c - even, s - odd)
        RealMatrix CMatrixEven = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());
        RealMatrix CMatrixOdd = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());

        for (int n = 0; n < psi.getNumber(); n++) {
            double lambdaEvenN = lambdasEven[n];
            double lambdaOddN = lambdasOdd[n];

            double AEvenN = AEven[n];
            double AOddN = AOdd[n];

            double lambdaEvenNSquare = lambdaEvenN * lambdaEvenN;
            double lambdaOddNSquare = lambdaOddN * lambdaOddN;

            double AEvenNSquare = AEvenN * AEvenN;
            double AOddNSquare = AOddN * AOddN;

            double lambdaEvenNCos = cos(lambdaEvenN);
            double lambdaOddNCos = cos(lambdaOddN);

            double lambdaEvenNSin = sin(lambdaEvenN);
            double lambdaOddNSin = sin(lambdaOddN);

            for (int m = 0; m < psi.getNumber(); m++) {
                double lambdaEvenM = lambdasEven[m];
                double lambdaOddM = lambdasOdd[m];

                double lambdaEvenMSquare = lambdaEvenM * lambdaEvenM;
                double lambdaOddMSquare = lambdaOddM * lambdaOddM;

                double AEvenM = AEven[m];
                double AOddM = AOdd[m];

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
                    CEven *= lambdaEvenM * lambdaEvenNCos * lambdaEvenMSin
                        - lambdaEvenN * lambdaOddMCos * lambdaEvenNSin;

                    COdd = 8 * lambdaOddNSquare * lambdaOddMSquare * AOddN * AOddM;
                    COdd /= lambdaOddNSquare * lambdaOddNSquare - lambdaOddMSquare * lambdaOddMSquare;
                    COdd *= lambdaOddN * lambdaOddNCos * lambdaOddMSin - lambdaOddM * lambdaEvenMCos * lambdaOddNSin;
                }

                CMatrixEven.addToEntry(n, m, CEven);
                CMatrixOdd.addToEntry(n, m, COdd);
            }
        }

        // 4. Calculate M matrices for even and odd cases (Mc[j], Ms[j], where c - even, s - odd)
        List<RealMatrix> MMatricesEven = new ArrayList<>();
        List<RealMatrix> MMatricesOdd = new ArrayList<>();

        int j = 10; // j - manually chosen parameter

        RealMatrix UJMatrixEven = new Array2DRowRealMatrix(j, ksi.getStepsNumber());
        RealMatrix UJMatrixOdd = new Array2DRowRealMatrix(j, ksi.getStepsNumber());

        RealMatrix PsiJMatrixEven = new Array2DRowRealMatrix(j, psi.getNumber());
        RealMatrix PsiJMatrixOdd = new Array2DRowRealMatrix(j, psi.getNumber());

        for (int jIndex = 0; jIndex < j; jIndex++) {
            for (int psiIndex = 0; psiIndex < psi.getNumber(); psiIndex++) {
                double lambdaEven = lambdasEven[psiIndex];
                double lambdaOdd = lambdasOdd[psiIndex];

                double AEvenValue = AEven[psiIndex];
                double AOddValue = AOdd[psiIndex];

                double PsiJNEven = pow(-1, jIndex) * 4 * pow(lambdaEven, 3) * AEvenValue * sin(lambdaEven)
                    / (pow(lambdaEven, 4) - pow(PI * jIndex, 4));

                PsiJMatrixEven.addToEntry(jIndex, psiIndex, PsiJNEven);

                if (0 == jIndex) {
                    continue;
                }

                double PsiJNOdd = pow(-1, jIndex) * (-4 * pow(lambdaOdd, 3)) * AOddValue * cos(lambdaOdd)
                    / (pow(lambdaOdd, 4) - pow(PI * jIndex - PI / 2, 4));

                PsiJMatrixOdd.addToEntry(jIndex, psiIndex, PsiJNOdd);
            }

            int ksiIndex = 0;
            for (double ksiValue : ksiValues) {
                double ksiSquared = ksiValue * ksiValue;

                double uEven = sqrt(ksiSquared + pow(PI * jIndex, 2));
                double uOdd = sqrt(ksiSquared + pow(PI * jIndex - PI / 2, 2));

                UJMatrixEven.addToEntry(jIndex, ksiIndex, uEven);
                UJMatrixOdd.addToEntry(jIndex, ksiIndex, uOdd);

                ksiIndex++;
            }
        }

        int ksiIndex = 0;
        for (double ksiValue : ksiValues) {
            RealMatrix MMatrixEven = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());
            RealMatrix MMatrixOdd = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());

            double ksiSquared = ksiValue * ksiValue;
            for (int m = 0; m < psi.getNumber(); m++) {
                for (int n = 0; n < psi.getNumber(); n++) {
                    double MEven = PsiJMatrixEven.getEntry(0, m) * PsiJMatrixEven.getEntry(0, n)
                        / (2 * ksiValue * tanh(ksiValue * physic.getHDimensionless()));

                    double MOdd = 0.0;

                    for (int jIndex = 1; jIndex < j; jIndex++) {
                        double ujEven = UJMatrixEven.getEntry(jIndex, ksiIndex);
                        double ujOdd = UJMatrixOdd.getEntry(jIndex, ksiIndex);

                        MEven += PsiJMatrixEven.getEntry(jIndex, m) * PsiJMatrixEven.getEntry(jIndex, n)
                            / ujEven * tanh(ujEven * physic.getHDimensionless());

                        MOdd += PsiJMatrixOdd.getEntry(jIndex, m) * PsiJMatrixOdd.getEntry(jIndex, n)
                            / ujOdd * tanh(ujOdd * physic.getHDimensionless());
                    }

                    MEven *= ksiSquared;
                    MOdd *= ksiSquared;

                    if (0.0 == ksiValue) {
                        MEven = (PsiJMatrixEven.getEntry(0, m) * PsiJMatrixEven.getEntry(0, n)) / (2 * physic.getHDimensionless());
                    }

                    MMatrixEven.addToEntry(m, n, MEven);
                    MMatrixOdd.addToEntry(m, n, MOdd);
                }
            }

            MMatricesEven.add(MMatrixEven);
            MMatricesOdd.add(MMatrixOdd);

            ksiIndex++;
        }

        // 5. Calculate D matrices for even and odd cases (Dc[j], Ds[j], where c - even, s - odd)
        List<RealMatrix> DMatricesEven = new ArrayList<>();
        List<RealMatrix> DMatricesOdd = new ArrayList<>();

        for (double ksiValue : ksiValues) {
            RealMatrix DMatrixEven = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());
            RealMatrix DMatrixOdd = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());

            double ksiFourth = pow(ksiValue, 4);
            for (int i = 0; i < psi.getNumber(); i++) {
                double DEven = pow(lambdasEven[i], 4) + ksiFourth;
                double DOdd = pow(lambdasOdd[i], 4) + ksiFourth;

                DMatrixEven.addToEntry(i, i, DEven);
                DMatrixOdd.addToEntry(i, i, DOdd);
            }

            DMatricesEven.add(DMatrixEven);
            DMatricesOdd.add(DMatrixOdd);
        }

        // 6. Calculate Q matrices for even and odd cases (Qc[j], Qs[j], where c - even, s - odd)
        List<RealMatrix> QMatricesEven = new ArrayList<>();
        List<RealMatrix> QMatricesOdd = new ArrayList<>();

        ksiIndex = 0;
        for (double ksiValue : ksiValues) {
            double ksiSquaredDouble = 2 * (ksiValue * ksiValue);

            RealMatrix DMatrixEven = DMatricesEven.get(ksiIndex);
            RealMatrix DMatrixOdd = DMatricesOdd.get(ksiIndex);

            RealMatrix KsiCMatrixEven = CMatrixEven.scalarMultiply(ksiSquaredDouble);
            RealMatrix KsiCMatrixOdd = CMatrixOdd.scalarMultiply(ksiSquaredDouble);

            RealMatrix QMatrixEven = DMatrixEven.subtract(KsiCMatrixEven);
            RealMatrix QMatrixOdd = DMatrixOdd.subtract(KsiCMatrixOdd);

            QMatricesEven.add(QMatrixEven);
            QMatricesOdd.add(QMatrixOdd);

            ksiIndex++;
        }

        // 7. Calculate R matrices for even and odd cases
        List<RealMatrix> RMatricesEven = new ArrayList<>();
        List<RealMatrix> RMatricesOdd = new ArrayList<>();

        ksiIndex = 0;
        for (double ksiValue : ksiValues) {
            RealMatrix QMatrixEven = QMatricesEven.get(ksiIndex);
            RealMatrix QMatrixOdd = QMatricesOdd.get(ksiIndex);

            RealMatrix MMatrixEven = MMatricesEven.get(ksiIndex);
            RealMatrix MMatrixOdd = MMatricesOdd.get(ksiIndex);

            double hFrSquared = physic.getHDimensionless() * physic.getFr() * physic.getFr();

            RealMatrix IMatrix = new Array2DRowRealMatrix(psi.getNumber(), psi.getNumber());
            for (int n = 0; n < psi.getNumber(); n++) {
                IMatrix.addToEntry(n, n, 1 - physic.getAlpha() * hFrSquared * ksiValue * ksiValue);
            }

            RealMatrix QMatrixEvenMultiplied = QMatrixEven.scalarMultiply(physic.getBeta());
            RealMatrix QMatrixOddMultiplied = QMatrixOdd.scalarMultiply(physic.getBeta());

            RealMatrix MMatrixEvenMultiplied = MMatrixEven.scalarMultiply(hFrSquared);
            RealMatrix MMatrixOddMultiplied = MMatrixOdd.scalarMultiply(hFrSquared);

            RealMatrix RMatrixEven = IMatrix.add(QMatrixEvenMultiplied).subtract(MMatrixEvenMultiplied);
            RealMatrix RMatrixOdd = IMatrix.add(QMatrixOddMultiplied).subtract(MMatrixOddMultiplied);

            RMatricesEven.add(RMatrixEven);
            RMatricesOdd.add(RMatrixOdd);

            ksiIndex++;
        }

        // 8. Calculate P1F arrays
        double[] P1FArray = new double[ksi.getStepsNumber()];

        double integrationLimitLower = 0.0;
        double integrationLimitUpper = 0.2;
        int stepsNumber = 320;

        ksiIndex = 0;
        for (double ksiValue : ksiValues) {
            TrapezedMethod trapezedMethod = new TrapezedMethod(
                (xn) -> P1.getYValue(xn) * cos(ksiValue * xn),
                integrationLimitLower,
                integrationLimitUpper,
                stepsNumber
            );

            P1FArray[ksiIndex] = sqrt(2 / PI) * trapezedMethod.solve();
            ksiIndex++;
        }

        // 9. Calculate Pm* array
        double[] PMArrayEven = new double[psi.getNumber()];
        double[] PMArrayOdd = new double[psi.getNumber()];

        for (int psiIndex = 0; psiIndex < psi.getNumber(); psiIndex++) {
            int psiIndexFinal = psiIndex;

            TrapezedMethod trapezedMethodEven = new TrapezedMethod(
                (xn) -> P2.getYValue(xn) * psiMatrixEven.get(psiIndexFinal).get(xn),
                P2.getXMin(),
                P2.getXMax(),
                P2.size() - 1
            );

            TrapezedMethod trapezedMethodOdd = new TrapezedMethod(
                (xn) -> P2.getYValue(xn) * psiMatrixOdd.get(psiIndexFinal).get(xn),
                P2.getXMin(),
                P2.getXMax(),
                P2.size() - 1
            );

            PMArrayEven[psiIndex] = trapezedMethodEven.solve();
            PMArrayOdd[psiIndex] = trapezedMethodOdd.solve();
        }

        // 10. Calculate P vectors for even and odd cases (Pm(ksi))
        List<RealMatrix> PVectorsEven = new ArrayList<>();
        List<RealMatrix> PVectorsOdd = new ArrayList<>();

        for (double P1F : P1FArray) {
            RealMatrix PVectorEven = new Array2DRowRealMatrix(psi.getNumber(), 1);
            RealMatrix PVectorOdd = new Array2DRowRealMatrix(psi.getNumber(), 1);

            for (int k = 0; k < psi.getNumber(); k++) {
                PVectorEven.addToEntry(k, 0, P1F * PMArrayEven[k]);
                PVectorOdd.addToEntry(k, 0, P1F * PMArrayOdd[k]);
            }

            PVectorsEven.add(PVectorEven);
            PVectorsOdd.add(PVectorOdd);
        }

        RealMatrix aRMatrixEven = new Array2DRowRealMatrix(ksi.getStepsNumber(), psi.getNumber());
        RealMatrix aRMatrixOdd = new Array2DRowRealMatrix(ksi.getStepsNumber(), psi.getNumber());
        RealMatrix aIMatrixEven = new Array2DRowRealMatrix(ksi.getStepsNumber(), psi.getNumber());
        RealMatrix aIMatrixOdd = new Array2DRowRealMatrix(ksi.getStepsNumber(), psi.getNumber());

        ksiIndex = 0;
        for (double ksiValue : ksiValues) {
            // 11.1 Calculate RFinal matrix for even and odd cases
            RealMatrix RMatrixEven = RMatricesEven.get(ksiIndex);
            RealMatrix RMatrixOdd = RMatricesOdd.get(ksiIndex);

            RealMatrix QMatrixEven = QMatricesEven.get(ksiIndex);
            RealMatrix QMatrixOdd = QMatricesOdd.get(ksiIndex);

            RealMatrix RMatrixEvenInverse = MatrixUtils.inverse(RMatrixEven);
            RealMatrix RMatrixOddInverse = MatrixUtils.inverse(RMatrixOdd);

            RealMatrix QRInverseMatrixEven = QMatrixEven.multiply(RMatrixEvenInverse);
            RealMatrix QRInverseMatrixOdd = QMatrixEven.multiply(RMatrixOddInverse);

            RealMatrix QSquaredRInverseMatrixEven = QRInverseMatrixEven.multiply(QMatrixEven);
            RealMatrix QSquaredRInverseMatrixOdd = QRInverseMatrixOdd.multiply(QMatrixOdd);

            double betaEpsilonKsi = physic.getBeta() * physic.getEpsilon() * ksiValue;
            double betaSquaredEpsilonKsiSquaredE = betaEpsilonKsi * physic.getBeta() * physic.getEpsilon() * ksiValue;

            RealMatrix RFinalMatrixEven = RMatrixEven.add(
                QSquaredRInverseMatrixEven.scalarMultiply(betaSquaredEpsilonKsiSquaredE)
            );

            RealMatrix RFinalMatrixOdd = RMatrixOdd.add(
                QSquaredRInverseMatrixOdd.scalarMultiply(betaSquaredEpsilonKsiSquaredE)
            );

            // 11.2 Calculate RFinal inverse negative matrix for even and odd cases
            RealMatrix RFinalMatrixInverseNegativeEven = MatrixUtils.inverse(RFinalMatrixEven).scalarMultiply(-1.0);
            RealMatrix RFinalMatrixInverseNegativeOdd = MatrixUtils.inverse(RFinalMatrixOdd).scalarMultiply(-1.0);

            // 11.3 Calculate ar for even and odd cases
            RealMatrix aRVectorEven = RFinalMatrixInverseNegativeEven.multiply(PVectorsEven.get(ksiIndex));
            RealMatrix aRVectorOdd = RFinalMatrixInverseNegativeOdd.multiply(PVectorsOdd.get(ksiIndex));

            // 11.4 Calculate al for even and odd cases
            RealMatrix RInverseQMatrixEven = RMatrixEvenInverse.multiply(QMatrixEven);
            RealMatrix RInverseQMatrixOdd = RMatrixOddInverse.multiply(QMatrixEven);

            RealMatrix aIVectorEven = RInverseQMatrixEven
                .multiply(aRVectorEven)
                .scalarMultiply(betaEpsilonKsi);

            RealMatrix aIVectorOdd = RInverseQMatrixOdd
                .multiply(aRVectorOdd)
                .scalarMultiply(betaEpsilonKsi);

            for (int i = 0; i < psi.getNumber(); i++) {
                aRMatrixEven.addToEntry(ksiIndex, i, aRVectorEven.getEntry(i, 0));
                aRMatrixOdd.addToEntry(ksiIndex, i, aRVectorOdd.getEntry(i, 0));
                aIMatrixEven.addToEntry(ksiIndex, i, aIVectorEven.getEntry(i, 0));
                aIMatrixOdd.addToEntry(ksiIndex, i, aIVectorOdd.getEntry(i, 0));
            }

            ksiIndex++;
        }

        // 12. Calculate W matrices for even and odd cases
        List<DeflectionPoint> points = new ArrayList<>();
        for (double xValue : x) {
            double[] wEvenArray = new double[psi.getNumber()];
            double[] wOddArray = new double[psi.getNumber()];

            for (int i = 0; i < psi.getNumber(); i++) {
                double QREven = 0.0d;
                double QIEven = 0.0d;

                double QROdd = 0.0d;
                double QIOdd = 0.0d;

                for (ksiIndex = 1; ksiIndex < ksi.getStepsNumber(); ksiIndex++) {
                    int ksiPreviousIndex = ksiIndex - 1;

                    double ksiPrevious = ksiValues[ksiPreviousIndex];
                    double ksiCurrent = ksiValues[ksiIndex];

                    double aRPreviousEven = aRMatrixEven.getEntry(ksiPreviousIndex, i);
                    double aIPreviousEven = aIMatrixEven.getEntry(ksiPreviousIndex, i);
                    double aRPreviousOdd = aRMatrixOdd.getEntry(ksiPreviousIndex, i);
                    double aIPreviousOdd = aIMatrixOdd.getEntry(ksiPreviousIndex, i);

                    double aRCurrentEven = aRMatrixEven.getEntry(ksiIndex, i);
                    double aICurrentEven = aIMatrixEven.getEntry(ksiIndex, i);
                    double aRCurrentOdd = aRMatrixOdd.getEntry(ksiIndex, i);
                    double aICurrentOdd = aIMatrixOdd.getEntry(ksiIndex, i);

                    double ARiEven = aRPreviousEven - ((aRCurrentEven - aRPreviousEven) / ksi.getStep()) * ksiPrevious;
                    double BRiEven = ((aRCurrentEven - aRPreviousEven) / ksi.getStep()) * ksiPrevious;

                    double ARiOdd = aRPreviousOdd - ((aRCurrentOdd - aRPreviousOdd) / ksi.getStep()) * ksiPrevious;
                    double BRiOdd = ((aRCurrentOdd - aRPreviousOdd) / ksi.getStep()) * ksiPrevious;

                    double AIiEven = aIPreviousEven - ((aICurrentEven - aIPreviousEven) / ksi.getStep()) * ksiPrevious;
                    double BIiEven = ((aICurrentEven - aIPreviousEven) / ksi.getStep()) * ksiPrevious;

                    double AIiOdd = aIPreviousOdd - ((aICurrentOdd - aIPreviousOdd) / ksi.getStep()) * ksiPrevious;
                    double BIiOdd = ((aICurrentOdd - aIPreviousOdd) / ksi.getStep()) * ksiPrevious;

                    double QRiEven;
                    double QRiOdd;

                    double QIiEven;
                    double QIiOdd;

                    if (xValue == 0.0d) {
                        QRiEven = ARiEven * ksi.getStep() + (BRiEven / 2.0d) * (ksiCurrent * ksiCurrent - ksiPrevious * ksiPrevious);
                        QRiOdd = ARiOdd * ksi.getStep() + (BRiOdd / 2.0d) * (ksiCurrent * ksiCurrent - ksiPrevious * ksiPrevious);

                        QIiEven = 0.0d;
                        QIiOdd = 0.0d;
                    } else {
                        QRiEven = (ARiEven + BRiEven * ksiCurrent) * (sin(ksiCurrent * xValue) / xValue)
                            - (ARiEven + BRiEven * ksiPrevious) * (sin(ksiPrevious * xValue) / xValue)
                            + BRiEven * (cos(ksiCurrent * xValue) / (xValue * xValue) - cos(ksiPrevious * xValue) / (xValue * xValue));

                        QRiOdd = (ARiOdd + BRiOdd * ksiCurrent) * (sin(ksiCurrent * xValue) / xValue)
                            - (ARiOdd + BRiOdd * ksiPrevious) * (sin(ksiPrevious * xValue) / xValue)
                            + BRiOdd * (cos(ksiCurrent * xValue) / (xValue * xValue) - cos(ksiPrevious * xValue) / (xValue * xValue));

                        QIiEven = (AIiEven + BIiEven * ksiPrevious) * (cos(ksiPrevious * xValue) / xValue)
                            - (AIiEven + BIiEven * ksiCurrent) * (cos(ksiCurrent * xValue) / xValue)
                            + BIiEven * (sin(ksiCurrent * xValue) / (xValue * xValue) - sin(ksiPrevious * xValue) / (xValue * xValue));

                        QIiOdd = (AIiOdd + BIiOdd * ksiPrevious) * (cos(ksiPrevious * xValue) / xValue)
                            - (AIiOdd + BIiOdd * ksiCurrent) * (cos(ksiCurrent * xValue) / xValue)
                            + BIiOdd * (sin(ksiCurrent * xValue) / (xValue * xValue) - sin(ksiPrevious * xValue) / (xValue * xValue));
                    }

                    QREven += QRiEven;
                    QIEven += QIiEven;

                    QROdd += QRiOdd;
                    QIOdd += QIiOdd;
                }

                wEvenArray[i] = QREven - QIEven;
                wOddArray[i] = QROdd - QIOdd;
            }

            for (double y : P2.getX()) {
                double z = 0.0d;
                for (int i = 0; i < psi.getNumber(); i++) {
                    double psiEven = psiMatrixEven.get(i).get(y);
                    double psiOdd = psiMatrixOdd.get(i).get(y);

                    z += psiEven * wEvenArray[i] + psiOdd * wOddArray[i];
                }

                z *= sqrt(2 / PI);

                points.add(new DeflectionPoint(xValue, y, z));
            }
        }

        return points;
    }

    /**
     * Calculates B even array
     *
     * @param lambdas a lambdas even array
     * @return a B even array
     */
    private double[] calculateBEven(double[] lambdas) {
        double[] B = new double[lambdas.length];
        for (int i = 0; i < lambdas.length; i++) {
            B[i] = cos(lambdas[i]) / cosh(lambdas[i]);
        }

        return B;
    }

    /**
     * Calculates B odd array
     *
     * @param lambdas a lambdas odd array
     * @return a B odd array
     */
    private double[] calculateBOdd(double[] lambdas) {
        double[] B = new double[lambdas.length];
        for (int i = 0; i < lambdas.length; i++) {
            B[i] = sin(lambdas[i]) / sinh(lambdas[i]);
        }

        return B;
    }

    /**
     * Calculates A even array
     *
     * @param B a B even array
     * @return an A even array
     */
    private double[] calculateAEven(double[] B) {
        double[] A = new double[B.length];
        for (int i = 0; i < B.length; i++) {
            A[i] = sqrt(1.0d / (1.0d + B[i] * B[i]));
        }

        return A;
    }

    /**
     * Calculates A odd array
     *
     * @param B a B odd array
     * @return an A odd array
     */
    private double[] calculateAOdd(double[] B) {
        double[] A = new double[B.length];
        for (int i = 0; i < B.length; i++) {
            A[i] = sqrt(1.0d / (1.0d - B[i] * B[i]));
        }

        return A;
    }
}
