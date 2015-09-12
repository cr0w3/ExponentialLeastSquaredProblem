import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;

public class Test {

	static int l;
	static double alpha, beta, gamma;
	static Double[] table;
	static double[] observedValues;

	public static void main(String[] args) {
		table = new Double[] { 362.0, 385.0, 432.0, 341.0, 382.0, 409.0, 498.0, 387.0, 473.0, 513.0, 582.0, 474.0,
				544.0, 582.0, 681.0, 557.0, 628.0, 707.0, 773.0, 592.0, 627.0, 725.0, 854.0, 661.0 };

		/*
		 * alpha = 0.7556; beta = 0.0; gamma = 0.9837;
		 */

		l = 4;

		optimizeParameters();

		// applyTripleExpSmoothing();

		// for (int i = 0; i < table.length; i++) {
		// System.out.println(table[i]);
		// }
	}

	private static void optimizeParameters() {

		final ExpParameterProblem problem = new ExpParameterProblem();
		problem.addTable(table);
		observedValues = new double[table.length - l];
		for (int i = l; i < table.length; i++) {
			problem.addPoint(table[i]);
			observedValues[i - l] = table[i];
		}

		final double[] start = { 0.5, 0.5, 0.5 };

		final Optimum optimum = getOptimizer()
				.optimize(builder(problem).target(observedValues).start(start).maxIterations(20).build());

		RealVector solution = optimum.getPoint();
		System.out.println("Solution: alpha = " + solution.getEntry(0) + ", beta = " + solution.getEntry(1) + ", gamma = " + solution.getEntry(2));
	}

	private static Double getInitialTrendTripleSmoothing() {
		if (table.length > 2 * l) {
			Double trend = 0D;
			for (int i = 0; i < l; i++) {
				trend += (Double.valueOf(table[i + l].toString()) - Double.valueOf(table[i].toString())) / l;
			}
			return trend / l;
		} else
			return -1D;
	}

	private static Double[] getSeasonalIndizes() {
		Double[] yearlyAverage = new Double[(int) Math.ceil(table.length / (double) l)];
		double sum = 0;
		for (int i = 1; i <= table.length; i++) {
			sum += Double.valueOf(table[i - 1].toString());
			if (i % l == 0) {
				yearlyAverage[(i / l) - 1] = sum / l;
				sum = 0;
			}
		}
		if (sum != 0)
			yearlyAverage[yearlyAverage.length - 1] = sum / (table.length % l);
		Double[] seasonalIndizes = new Double[l];
		double seasonalIndex = 0;
		for (int i = 0; i < l; i++) {
			for (int j = i; j < table.length; j += l) {
				seasonalIndex += Double.valueOf(table[j].toString()) / yearlyAverage[j / l];
			}
			int seasonalEntries = table.length / l;
			seasonalIndizes[i] = seasonalIndex / seasonalEntries;
			seasonalIndex = 0;
		}
		return seasonalIndizes;
	}

	private static void applyTripleExpSmoothing() {
		Double tripleWeightedAverage = Double.valueOf(table[l - 1].toString());
		Double tripleWeightedAveragePrev = 0D;
		Double trend = getInitialTrendTripleSmoothing();
		Double[] seasonalIndizes = getSeasonalIndizes();
		for (int i = l; i < table.length; i++) {
			tripleWeightedAveragePrev = tripleWeightedAverage;

			tripleWeightedAverage = alpha * Double.valueOf(table[i].toString()) / seasonalIndizes[i % l]
					+ (1 - alpha) * (tripleWeightedAveragePrev + trend);

			trend = beta * (tripleWeightedAverage - tripleWeightedAveragePrev) + (1 - beta) * trend;

			double index = seasonalIndizes[i % l];
			seasonalIndizes[i % l] = gamma * Double.valueOf(table[i].toString()) / tripleWeightedAverage
					+ (1 - gamma) * seasonalIndizes[i % l];

			table[i] = (tripleWeightedAverage + trend) * index;
		}
	}

	private static class ExpParameterProblem {

		private List<Double> count;
		private Double[] table;
		private HelperClass helper;

		public ExpParameterProblem() {
			count = new ArrayList<Double>();
		}

		public void addTable(Double[] table) {
			this.table = table;
		}

		public void addPoint(Double p) {
			count.add(p);
		}

		public MultivariateVectorFunction getModelFunction() {
			return new MultivariateVectorFunction() {
				public double[] value(double[] params) {
					helper = new HelperClass(table, l);
					double[] values = new double[count.size()];
					for (int i = 0; i < values.length; ++i) {
						helper.calculateSmoothedObservation(i, params[0], params[1], params[2]);
						double modelF = (helper.getSmoothedObservation() + helper.getTrend())
								* helper.getPrevSeasonalIndizes()[i % l];
						values[i] = modelF;
					}
					return values;
				}
			};
		}

		public MultivariateMatrixFunction getModelFunctionJacobian() {
			return new MultivariateMatrixFunction() {
				public double[][] value(double[] params) {
					double[][] jacobian = new double[observedValues.length][3];
					helper = new HelperClass(table, l);
					for (int i = 0; i < jacobian.length; ++i) {
						double prevPrevSeasonalIndex = helper.getPrevSeasonalIndizes()[i % l];
						helper.calculateSmoothedObservation(i, params[0], params[1], params[2]);
						jacobian[i][0] = table[i + l] / helper.getPrevSeasonalIndizes()[i % l]
								- helper.getPrevTrend() - helper.getPrevSmoothedObservation();
						jacobian[i][1] = helper.getSmoothedObservation() - helper.getPrevSmoothedObservation()
								- helper.getPrevTrend();
						jacobian[i][2] = table[i] / helper.getPrevLSmoothedObservations()[i % l]
								- prevPrevSeasonalIndex;
					}
					return jacobian;
				}
			};
		}
	}

	public static LeastSquaresBuilder base() {
		return new LeastSquaresBuilder().checkerPair(new SimpleVectorValueChecker(1e-6, 1e-6)).maxEvaluations(1000)
				.maxIterations(getMaxIterations());
	}

	public static LeastSquaresBuilder builder(ExpParameterProblem problem) {
		return base().model(problem.getModelFunction(), problem.getModelFunctionJacobian());
	}

	public static int getMaxIterations() {
		return 1000;
	}

	public static LeastSquaresOptimizer getOptimizer() {
		return new LevenbergMarquardtOptimizer();
	}
}
