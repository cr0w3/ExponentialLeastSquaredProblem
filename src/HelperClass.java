
public class HelperClass {

	private Double[] observations;

	private Double[] seasonalIndizes;
	private Double[] prevSeasonalIndizes;

	private Double prevTrend;
	private Double prevSmoothedObservation;
	private Double[] prevLSmoothedObservations;

	private Double smoothedObservation;
	private Double trend;
	
	private int l;

	public HelperClass(Double[] table, int l) {
		this.l = l;
		seasonalIndizes = new Double[l];
		prevSeasonalIndizes = new Double[l];
		prevLSmoothedObservations = new Double[l];
		this.observations = table;
		for (int i = 0; i < l; i++) {
			prevLSmoothedObservations[i] = observations[i];
		}
		smoothedObservation = observations[l - 1];
		prevSmoothedObservation = smoothedObservation;
		calcInitialTrendTripleSmoothing();
		calcSeasonalIndizes();
	}

	private void calcInitialTrendTripleSmoothing() {
		if (observations.length >= 2 * l) {
			Double trend = 0D;
			for (int i = 0; i < l; i++) {
				trend += (Double.valueOf(observations[i + l].toString()) - Double.valueOf(observations[i].toString()))
						/ l;
			}
			prevTrend = trend / l;
			this.trend = prevTrend;
		}
	}

	private void calcSeasonalIndizes() {
		Double[] yearlyAverage = new Double[(int) Math.ceil(observations.length / (double) l)];
		double sum = 0;
		for (int i = 1; i <= observations.length; i++) {
			sum += Double.valueOf(observations[i - 1].toString());
			if (i % l == 0) {
				yearlyAverage[(i / l) - 1] = sum / l;
				sum = 0;
			}
		}
		if (sum != 0)
			yearlyAverage[yearlyAverage.length - 1] = sum / (observations.length % l);
		Double[] seasonalIndizes = new Double[l];
		double seasonalIndex = 0;
		for (int i = 0; i < l; i++) {
			for (int j = i; j < observations.length; j += l) {
				seasonalIndex += Double.valueOf(observations[j].toString()) / yearlyAverage[j / l];
			}
			int seasonalEntries = observations.length / l;
			seasonalIndizes[i] = seasonalIndex / seasonalEntries;
			seasonalIndex = 0;
		}
		this.seasonalIndizes = seasonalIndizes;
		this.prevSeasonalIndizes = seasonalIndizes;
	}

	public void calculateSmoothedObservation(int index, double a, double b, double c) {
		if (index >= l) {
			prevLSmoothedObservations[index % l] = smoothedObservation;
		}
		prevSmoothedObservation = smoothedObservation;
		smoothedObservation = a * observations[index + l] / seasonalIndizes[index % l]
				+ (1 - a) * (prevSmoothedObservation + prevTrend);
		prevTrend = trend;
		trend = b * (smoothedObservation - prevSmoothedObservation) + (1 - b) * prevTrend;
		prevSeasonalIndizes[index % l] = seasonalIndizes[index % l];
		seasonalIndizes[index % l] = c * observations[index + l] / smoothedObservation
				+ (1 - c) * seasonalIndizes[index % l];
	}

	public Double getPrevTrend() {
		return prevTrend;
	}

	public Double[] getPrevSeasonalIndizes() {
		return prevSeasonalIndizes;
	}

	public Double[] getSeasonalIndizes() {
		return seasonalIndizes;
	}

	public Double getTrend() {
		return trend;
	}

	public Double getPrevSmoothedObservation() {
		return prevSmoothedObservation;
	}

	public Double getSmoothedObservation() {
		return smoothedObservation;
	}

	public Double[] getPrevLSmoothedObservations() {
		return prevLSmoothedObservations;
	}
}
