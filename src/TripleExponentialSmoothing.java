
/**
 * See <a href=
 * 'http://www.itl.nist.gov/div898/handbook/pmc/section4/pmc435.htm'>this
 * link</a> for a better understanding of triple exponential smoothing.
 * 
 * <p>Equations: <ul> <li>S(t) = alpha * X(t) / I (t-L) + (1 - alpha)(S(t-1) +
 * b(t-1))</li> <li>b(t) = beta * (S(t) - S(t-1)) + (1 - beta) * b(t-1)</li>
 * <li>I(t) = gamma * X(t) / S(t) + (1 - gamma) * I(t-L)</li> <li>F(t) = (S(t) +
 * b(t)) * I(t-L)</li></ul> Please notice the usage of the variables equals the
 * following parameters: <ul> <li>S(t) = {@link #smoothedObservation}</li>
 * <li>S(t-1) = {@link #prevSmoothedObservation}</li> <li>b(t) = {@link #trend}
 * </li> <li>b(t-1) = {@link #prevTrend}</li> <li>I(t-L) =
 * {@link #seasonalIndizes} </li> <li>I(t-2L) = {@link #prevSeasonalIndizes}</li>
 * <li>F(t) = The result of
 * {@link #calculateSmoothedObservation(int, double, double, double)}</li> </ul>
 * </p>
 * 
 * @author cr0w3
 *
 */
public class TripleExponentialSmoothing {

	/**
	 * The actual observations X(0) ... X(t)
	 */
	private double[] observations;

	/**
	 * The seasonal indizes of length valuesPerSeason. We always work with
	 * I(t-L) so we store them until we reach the next cycle of seasons.
	 */
	private double[] seasonalIndizes;
	/**
	 * A storage for previously used seasonal indizes, i.e., I(t-2L)
	 */
	private double[] prevSeasonalIndizes;

	/**
	 * The previously used trend b(t-1)
	 */
	private double prevTrend;
	/**
	 * The previous smoothed observation S(t-1)
	 */
	private double prevSmoothedObservation;
	/**
	 * A storage of previous smoothed observations used for the derivation of
	 * I(t-L)
	 */
	private double[] prevLSmoothedObservations;

	/**
	 * The current smoothed observation S(t)
	 */
	private double smoothedObservation;
	/**
	 * The current trend b(t)
	 */
	private double trend;

	/**
	 * The L value in the formula
	 */
	private int valuesPerSeason;

	/**
	 * Constructs a new {@link TripleExponentialSmoothing} Object.
	 * 
	 * @param table
	 *            The complete data table including all observations
	 * @param valuesPerSeason
	 *            The L value, i.e., how many values are recorded per season
	 */
	public TripleExponentialSmoothing(double[] table, int valuesPerSeason) {
		this.valuesPerSeason = valuesPerSeason;
		seasonalIndizes = new double[valuesPerSeason];
		prevSeasonalIndizes = new double[valuesPerSeason];
		prevLSmoothedObservations = new double[valuesPerSeason];
		this.observations = table;
		// we fill prevLSmoothedObservations with the first L values
		// because we do not calculate an S(t) for them.
		for (int i = 0; i < valuesPerSeason; i++) {
			prevLSmoothedObservations[i] = observations[i];
		}
		// the first S(t-1) is S(L-1)
		smoothedObservation = observations[valuesPerSeason - 1];
		prevSmoothedObservation = smoothedObservation;
		calcInitialTrendTripleSmoothing();
		calcSeasonalIndizes();
	}

	/**
	 * Calculate the initial trend of triple exp smoothing. Equation: 1 / L *
	 * ((X(1+L) - X(1))/L + (X(2+L) - X(2))/L + ... + (X(L+L) - X(L)) / L
	 */
	private void calcInitialTrendTripleSmoothing() {
		if (observations.length >= 2 * valuesPerSeason) {
			Double trend = 0D;
			for (int i = 0; i < valuesPerSeason; i++) {
				// Calculation of each (X(i+L) - X(i)) / L
				trend += (observations[i + valuesPerSeason] - observations[i]) / valuesPerSeason;
			}
			// The 1/L * (...) part.
			prevTrend = trend / valuesPerSeason;
			this.trend = prevTrend;
		}
	}

	/**
	 * Calculate the initial seasonal indizes. This consists of three steps:
	 * <ol>
	 * <li>Compute the average of each season: Average(season) = SUM(X(i)) /
	 * length of the season</li>
	 * <li>Form the seasonal indizes (consider L = 4): I(1) = (X(1)/Average(1) +
	 * X(5)/Average(2) + X(9)/Average(3) ...) / Number of Averages</li>
	 * </ol>
	 */
	private void calcSeasonalIndizes() {
		// The length of the season is calculated by dividing the observations
		// length / L
		double[] seasonalAverage = new double[(int) Math.ceil(observations.length / (double) valuesPerSeason)];
		double sum = 0;
		for (int i = 1; i <= observations.length; i++) {
			// We calculate the mean of each season
			sum += observations[i - 1];
			// Once L-values have been summed up
			if (i % valuesPerSeason == 0) {
				// we divide them through L to get the mean
				seasonalAverage[(i / valuesPerSeason) - 1] = sum / valuesPerSeason;
				// and start anew
				sum = 0;
			}
		}
		// If there are leftovers (i.e., 30 values but L = 4, we have 2 values
		// left that were not summed up
		if (sum != 0)
			// We use the sum of these 2 values and devide it by 2 (30 % 4) to
			// get the mean of the last season
			seasonalAverage[seasonalAverage.length - 1] = sum / (observations.length % valuesPerSeason);
		// We have L seasonal indizes to calculate
		double[] seasonalIndizes = new double[valuesPerSeason];
		double seasonalIndex = 0;
		// for each season
		for (int i = 0; i < valuesPerSeason; i++) {
			for (int j = i; j < observations.length; j += valuesPerSeason) {
				// We take the current observation and divide it by the seasonal
				// average and sum them up
				seasonalIndex += observations[j] / seasonalAverage[j / valuesPerSeason];
			}
			// Afterwards we divide the seasonal index through the number of
			// seasons
			seasonalIndizes[i] = seasonalIndex / seasonalAverage.length;
			seasonalIndex = 0;
		}
		this.seasonalIndizes = seasonalIndizes;
		this.prevSeasonalIndizes = seasonalIndizes;
	}

	/**
	 * The calculation of each equation. Index t = 0 will be treated as t = L as
	 * we do not calculate the smoothed observation for the first L values.
	 * 
	 * @param index
	 *            The current index (t) of X
	 * @param a
	 *            Alpha
	 * @param b
	 *            Beta
	 * @param c
	 *            Gamma
	 */
	public void calculateSmoothedObservation(int index, double a, double b, double c) {
		// if index is greater or equal than L
		if (index >= valuesPerSeason) {
			// we need to store prevLSmoothedObservations for the deviation
			prevLSmoothedObservations[index % valuesPerSeason] = smoothedObservation;
		}
		// S(t-1) = S(t) before
		prevSmoothedObservation = smoothedObservation;
		// calculating S(t) = alpha * X(t) / I(t-L) + (1 - alpha) * (S(t-1) + b(t-1)
		smoothedObservation = a * observations[index + valuesPerSeason] / seasonalIndizes[index % valuesPerSeason]
				+ (1 - a) * (prevSmoothedObservation + prevTrend);
		// b(t-1) = b(t) before
		prevTrend = trend;
		// calculating b(t) = beta * (S(t) - S(t-1)) + (1 - beta) * b(t-1)
		trend = b * (smoothedObservation - prevSmoothedObservation) + (1 - b) * prevTrend;
		// I(t-2L) = I(t-L)
		prevSeasonalIndizes[index % valuesPerSeason] = seasonalIndizes[index % valuesPerSeason];
		// before calculating I(t) = gamma * X(t) / S(t) + (1 - gamma) * I(t-L)
		seasonalIndizes[index % valuesPerSeason] = c * observations[index + valuesPerSeason] / smoothedObservation
				+ (1 - c) * seasonalIndizes[index % valuesPerSeason];
	}

	/**
	 * Returns b(t-1).
	 * @return b(t-1)
	 */
	public double getPrevTrend() {
		return prevTrend;
	}

	/**
	 * Returns I(t-2L).
	 * @return I(t-2L)
	 */
	public double[] getPrevSeasonalIndizes() {
		return prevSeasonalIndizes;
	}

	/**
	 * Returns I(t-L).
	 * @return I(t-L)
	 */
	public double[] getSeasonalIndizes() {
		return seasonalIndizes;
	}

	/**
	 * Returns b(t).
	 * @return b(t)
	 */
	public double getTrend() {
		return trend;
	}

	/**
	 * Returns S(t-1).
	 * @return S(t-1)
	 */
	public double getPrevSmoothedObservation() {
		return prevSmoothedObservation;
	}

	/**
	 * Returns S(t).
	 * @return S(t)
	 */
	public double getSmoothedObservation() {
		return smoothedObservation;
	}

	/**
	 * Returns S(t-L).
	 * @return S(t-L)
	 */
	public double[] getPrevLSmoothedObservations() {
		return prevLSmoothedObservations;
	}
}
