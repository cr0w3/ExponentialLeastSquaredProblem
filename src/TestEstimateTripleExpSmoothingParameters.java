import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

public class TestEstimateTripleExpSmoothingParameters {

	/**
	 * L value
	 */
	static int entriesPerSeason;
	static double alpha, beta, gamma;
	/**
	 * All observations
	 */
	static double[] table;
	/**
	 * All observations starting with X(L)
	 */
	static double[] observedValues;
	/**
	 * The helper that supports us with calculating all necessary values
	 */
	static TripleExponentialSmoothing helper;

	public static void main(String[] args) {
		// The example values of <a
		// href='http://www.itl.nist.gov/div898/handbook/pmc/section4/pmc436.htm'>this
		// page</a>
		table = new double[] { 362.0, 385.0, 432.0, 341.0, 382.0, 409.0, 498.0, 387.0, 473.0, 513.0, 582.0, 474.0,
				544.0, 582.0, 681.0, 557.0, 628.0, 707.0, 773.0, 592.0, 627.0, 725.0, 854.0, 661.0 };

		/*
		 * The correct parameters, used to test applyTripleExpSmoothing(); Not
		 * important here. alpha = 0.7556; beta = 0.0; gamma = 0.9837;
		 */

		// L = 4
		entriesPerSeason = 4;

		// try to optimize alpha, beta and gamma
		optimizeParameters();

		// applyTripleExpSmoothing();

		// for (int i = 0; i < table.length; i++) {
		// System.out.println(table[i]);
		// }
	}

	private static void optimizeParameters() {
		// we fill the observedValues array with all X values from X(L) to
		// X(24).
		observedValues = new double[table.length - entriesPerSeason];
		for (int i = entriesPerSeason; i < table.length; i++) {
			observedValues[i - entriesPerSeason] = table[i];
		}
		// we use a MultivariateJacobianFunction to calculate the distance
		// between F(t) and X(t)
		MultivariateJacobianFunction distancesToCurrentFValue = new MultivariateJacobianFunction() {
			public Pair<RealVector, RealMatrix> value(final RealVector point) {
				// retrieve our three parameters
				double alpha = point.getEntry(0);
				double beta = point.getEntry(1);
				double gamma = point.getEntry(2);

				// a vector with all observations (X(L), ...)
				RealVector value = new ArrayRealVector(observedValues.length);
				// a jacobian matrix with the deviations of F(t)
				RealMatrix jacobian = new Array2DRowRealMatrix(observedValues.length, 3);
				// instantiate a new helper with all observations and L
				helper = new TripleExponentialSmoothing(table, entriesPerSeason);
				for (int i = 0; i < observedValues.length; ++i) {
					// Set up our jacobian matrix:
					// F(t) = ((alpha * X(t) / I(t-L) + (1 - alpha) * (S(t-1) +
					// b(t-1)) + ... (first part)
					// derivative with respect to dF(t)/dalpha =
					// X(t) / I(t-L) - S(t-1) - b(t-1)
					jacobian.setEntry(i, 0, table[i] / helper.getSeasonalIndizes()[i % entriesPerSeason]
							- helper.getPrevSmoothedObservation() - helper.getPrevTrend());
					// F(t) = ... + beta * (S(t) - S(t-1) + (1-beta) * b(t-1)) *
					// ... (second part)
					// derivative with respect to dF(t)/dbeta =
					// S(t) - S(t-1) - b(t-1)
					jacobian.setEntry(i, 1, helper.getSmoothedObservation() - helper.getPrevSmoothedObservation()
							- helper.getPrevTrend());
					// F(t) = ... * I(t-L)
					// I(t) = gamma * X(t) / S(t) + (1 - gamma) * I(t-L)
					// thus, I(t-L) = gamma * X(t-L) / S(t-L) + (1 - gamma) *
					// I(t-2L)
					jacobian.setEntry(i, 2, table[i] / helper.getPrevLSmoothedObservations()[i % entriesPerSeason]
							- helper.getPrevSeasonalIndizes()[i % entriesPerSeason]);
					// update our TripleExponentialSmoothing values as we now
					// face the first entry after calculating our initial values
					helper.calculateSmoothedObservation(i, alpha, beta, gamma);
					// F(t) = (S(t) + b(t)) * I(t-L)
					double fValue = (helper.getSmoothedObservation() + helper.getTrend())
							* helper.getSeasonalIndizes()[i % entriesPerSeason];
					// set F(t)
					value.setEntry(i, fValue);
				}

				return new Pair<RealVector, RealMatrix>(value, jacobian);

			}
		};

		// we start with alpha = 0.2, beta = 0.3, gamma = 0.4
		final double[] start = { 0.2, 0.3, 0.4 };

		// build our problem
		LeastSquaresProblem problem = new LeastSquaresBuilder().start(start).model(distancesToCurrentFValue)
				.target(observedValues).lazyEvaluation(false).maxEvaluations(1000).maxIterations(1000).build();

		// optimize
		LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);

		// retrieve our solution
		RealVector solution = optimum.getPoint();
		
		// Gets us 0.2, 0.3 and 0.4, which does not work properly.
		System.out.println("Solution: alpha = " + solution.getEntry(0) + ", beta = " + solution.getEntry(1)
				+ ", gamma = " + solution.getEntry(2));
	}

	// Ignore these methods:
	/*private static double getInitialTrendTripleSmoothing() {
		if (table.length > 2 * entriesPerSeason) {
			double trend = 0D;
			for (int i = 0; i < entriesPerSeason; i++) {
				trend += (table[i + entriesPerSeason] - table[i]) / entriesPerSeason;
			}
			return trend / entriesPerSeason;
		} else
			return -1D;
	}

	// Ignore
	private static double[] getSeasonalIndizes() {
		double[] yearlyAverage = new double[(int) Math.ceil(table.length / (double) entriesPerSeason)];
		double sum = 0;
		for (int i = 1; i <= table.length; i++) {
			sum += table[i - 1];
			if (i % entriesPerSeason == 0) {
				yearlyAverage[(i / entriesPerSeason) - 1] = sum / entriesPerSeason;
				sum = 0;
			}
		}
		if (sum != 0)
			yearlyAverage[yearlyAverage.length - 1] = sum / (table.length % entriesPerSeason);
		double[] seasonalIndizes = new double[entriesPerSeason];
		double seasonalIndex = 0;
		for (int i = 0; i < entriesPerSeason; i++) {
			for (int j = i; j < table.length; j += entriesPerSeason) {
				seasonalIndex += table[j] / yearlyAverage[j / entriesPerSeason];
			}
			int seasonalEntries = table.length / entriesPerSeason;
			seasonalIndizes[i] = seasonalIndex / seasonalEntries;
			seasonalIndex = 0;
		}
		return seasonalIndizes;
	}

	// Ignore
	private static void applyTripleExpSmoothing() {
		double tripleWeightedAverage = table[entriesPerSeason - 1];
		double tripleWeightedAveragePrev = 0D;
		double trend = getInitialTrendTripleSmoothing();
		double[] seasonalIndizes = getSeasonalIndizes();
		for (int i = entriesPerSeason; i < table.length; i++) {
			tripleWeightedAveragePrev = tripleWeightedAverage;

			tripleWeightedAverage = alpha * table[i] / seasonalIndizes[i % entriesPerSeason]
					+ (1 - alpha) * (tripleWeightedAveragePrev + trend);

			trend = beta * (tripleWeightedAverage - tripleWeightedAveragePrev) + (1 - beta) * trend;

			double index = seasonalIndizes[i % entriesPerSeason];
			seasonalIndizes[i % entriesPerSeason] = gamma * table[i] / tripleWeightedAverage
					+ (1 - gamma) * seasonalIndizes[i % entriesPerSeason];

			table[i] = (tripleWeightedAverage + trend) * index;
		}
	} */
}
