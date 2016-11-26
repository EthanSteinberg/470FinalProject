#pragma once

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>

#include "mathutils.h"

/**
 * A class for holding a bunch of histories mapped to probabilities.
 */
class densemap {

public:

	/**
	 * Dummy constructor. Doesn't actually initialize it.
	 */
	densemap() {
		initialized = false;
	}

	/**
	 * Check if the map is initialized.
	 */
	bool isInitialized() const {
		return initialized;
	}

	/**
	 * Initialize the map.
	 * taxa_bits are the taxas in this map.
	 * netNodeChoices are the choices at each network node.
	 */
	void init(uint16_t taxa_bits, std::vector<int64_t> netNodeChoices) {
		initialized = true;
		std::memset(histories, 0, sizeof(histories));
		history_bitset = 0;
		this->taxa_bits = taxa_bits;
		choices = netNodeChoices;
	}

	/**
	 * Check if the two densmaps are compatabile (have the same choices)
	 */
	bool isCompatible(const densemap& other) const {
		for (unsigned int i = 0; i < choices.size(); i++) {
			if (choices[i] != other.choices[i] && choices[i] != -1 && other.choices[i] != -1) {
				return false;
			}
		}

		return true;
	}

	/**
	 * Add a value to the history.
	 */
	void addToHistory(int history, double value) {
		histories[history] += value;
		history_bitset |= 1LL << history;
	}

	/**
	 * Set a history value.
	 */
	void setHistory(int history, double value) {
		histories[history] = value;
		history_bitset |= 1LL << history;
	}

	/**
	 * Get a history value.
	 */
	double getHistory(int history) const {
		return histories[history];
	}

	/**
	 * Get a bitset for the histories available.
	 */
	uint64_t getHistoryBitset() const {
		return history_bitset;
	}

	/**
	 * Return the maximum amount of histories for a given number of events.
	 */
	uint64_t getMaxHistory(int numEvents) const {
		return 1 << numEvents;
	}

	/**
	 * Get the current taxa bits.
	 */
	uint16_t getTaxaBits() const {
		return taxa_bits;
	}

	/**
	 * Add two densmaps together.
	 */
	densemap& operator+=(const densemap& rhs) {
		uint64_t rhsBitset = rhs.getHistoryBitset();

		while (rhsBitset != 0) {
			int rhsOne = 63 - __builtin_clzll(rhsBitset);
			rhsBitset ^= (1LL << rhsOne);

			addToHistory(rhsOne, rhs.getHistory(rhsOne));
		}

		return *this;
	}


	/**
	 * The current choices.
	 */
	std::vector<int64_t> choices;

private:
	// If this map is initialized.
	bool initialized;

	// The current taxa bits.
	uint16_t taxa_bits;

	// All the histories.
	double histories[1 << 6];

	// Which histories are set.
	uint64_t history_bitset;

};

/**
 * Merge the choices of two densemaps.
 */
inline std::vector<int64_t> mergeChoices(const densemap& left, const densemap& right) {
	std::vector<int64_t> result(left.choices.size());

	for (unsigned int i = 0; i <result.size();i++) {
		if (left.choices[i] == -1) {
			result[i] = right.choices[i];
		} else if (right.choices[i] == -1) {
			result[i] = left.choices[i];
		} else {
			result[i] = left.choices[i];
		}
	}

	return result;
}

/**
 * Combine two densemaps.
 */
inline densemap combine(const densemap& left, const densemap& right) {
	densemap result;
	result.init(left.getTaxaBits() | right.getTaxaBits(), mergeChoices(left, right));

	uint64_t leftBitset = left.getHistoryBitset();

	while (leftBitset != 0) {
		int leftOne = 63 - __builtin_clzll(leftBitset);
		leftBitset ^= (1LL << leftOne);

		uint64_t rightBitset = right.getHistoryBitset();

		while (rightBitset != 0) {
			int rightOne = 63 - __builtin_clzll(rightBitset);
			rightBitset ^= (1LL << rightOne);

			result.setHistory(leftOne | rightOne, left.getHistory(leftOne) * right.getHistory(rightOne));
		}
	}

	return result;
}

/**
 * Combine the derivatives of densemaps.
 */
inline densemap combineDerivatives(const densemap& left, const densemap& leftDerivative, const densemap& right, const densemap& rightDerivative) {
	densemap result;
	result.init(left.getTaxaBits() | right.getTaxaBits(), mergeChoices(left, right));

	uint64_t leftBitset = left.getHistoryBitset();

	while (leftBitset != 0) {
		int leftOne = 63 - __builtin_clzll(leftBitset);
		leftBitset ^= (1LL << leftOne);

		uint64_t rightBitset = right.getHistoryBitset();

		while (rightBitset != 0) {
			int rightOne = 63 - __builtin_clzll(rightBitset);
			rightBitset ^= (1LL << rightOne);

			result.setHistory(leftOne | rightOne, leftDerivative.getHistory(leftOne) * right.getHistory(rightOne) + left.getHistory(leftOne) * rightDerivative.getHistory(rightOne));
		}
	}

	return result;
}

/**
 * Combine a list of densemaps.
 */
inline std::vector<densemap> combine(const std::vector<densemap>& left, const std::vector<densemap>& right) {
	std::vector<densemap> result;
	for (auto&& leftOne : left) {
		for (auto&& rightOne : right) {
			if (leftOne.isCompatible(rightOne)) {
				// You can only merge when the two share no taxa bits
				result.push_back(combine(leftOne, rightOne));
			}
		}
	}
	return result;
}

/**
 * Combine the derivatives for a list of densemaps.
 */
inline std::vector<densemap> combineDerivatives(const std::vector<densemap>& left, const std::vector<densemap>& leftDerivatives, const std::vector<densemap>& right, const std::vector<densemap>& rightDerivatives) {
	std::vector<densemap> result;
	for (unsigned int leftIndex = 0; leftIndex < left.size(); leftIndex ++) {
		for (unsigned int rightIndex = 0; rightIndex < right.size(); rightIndex ++) {
			auto&& leftOne = left[leftIndex];
			auto&& rightOne = right[rightIndex];
			if (leftOne.isCompatible(rightOne)) {
				// You can only merge when the two share no taxa bits
				result.push_back(combineDerivatives(leftOne, leftDerivatives[leftIndex], rightOne, rightDerivatives[rightIndex]));
			}
		}
	}
	return result;
}

/**
 * Perform BFS on a given history.
 */
inline std::pair<std::array<int, 1<<6 >, uint64_t> performBFS(uint8_t history, uint16_t taxaBits, const std::vector<int>& events) {
	std::array<int, 1<<6> numberOfWaysToReach = {};
	uint64_t numberOfWaysBitset = 0;

	int queue[1 << 6];
	int currentIndex = 0;

	queue[currentIndex++] = history;

	numberOfWaysToReach[history] = 1;
	numberOfWaysBitset |= 1LL << history;

	for (int i = 0; i < currentIndex; i++) {
		int next = queue[i];
		auto full = taxaBits | next;

		for (unsigned int i = 0; i < events.size(); i++) {
			auto& event = events[i];

			if ((full & event) == event && (next & 1 << i) == 0) {
				int final = next | 1 << i;

				if (numberOfWaysToReach[final] == 0) {
					queue[currentIndex++] = final;
				}
				numberOfWaysToReach[final] += numberOfWaysToReach[next];
				numberOfWaysBitset |= 1LL << final;
			}
		}
	}

	return {numberOfWaysToReach, numberOfWaysBitset };
}

/**
 * Update a densemap along a certain amount of time.
 */
inline densemap update(const densemap& current, const std::vector<int>& events, double length) {
	densemap result;
	result.init(current.getTaxaBits(), current.choices);

	uint64_t bitset = current.getHistoryBitset();
	while (bitset != 0) {
		int history = 63 - __builtin_clzll(bitset);
		bitset ^= (1LL << history);

		std::array<int, 1<<6> numberOfWaysToReach = {};
		uint64_t numberOfWaysBitset = 0;

		std::tie(numberOfWaysToReach, numberOfWaysBitset) = performBFS(history, current.getTaxaBits(), events);

		while (numberOfWaysBitset != 0) {
			int reachable = 63 - __builtin_clzll(numberOfWaysBitset);
			numberOfWaysBitset ^= (1LL << reachable);

			int numberOfWays = numberOfWaysToReach[reachable];

			int startingCount = __builtin_popcount(result.getTaxaBits()) - __builtin_popcount(history);
			int finalCount = __builtin_popcount(result.getTaxaBits()) - __builtin_popcount(reachable);

			double numberOfOptions = getNumberOfOptions(startingCount, finalCount);

			double weight = (double) numberOfWays / numberOfOptions;

			double total = current.getHistory(history) * weight * puv(startingCount, finalCount, length);

			if (total != 0) {
				result.addToHistory(reachable, total);
			}
		}

	}

	return result;
}

/**
 * Update a the derivative of a densemap along a certain amount of time.
 */
inline densemap derivativeUpdate(const densemap& current, const std::vector<int>& events, double length) {
	densemap result;
	result.init(current.getTaxaBits(), current.choices);

	uint64_t bitset = current.getHistoryBitset();
	while (bitset != 0) {
		int history = 63 - __builtin_clzll(bitset);
		bitset ^= (1LL << history);

		std::array<int, 1<<6> numberOfWaysToReach = {};
		uint64_t numberOfWaysBitset = 0;

		std::tie(numberOfWaysToReach, numberOfWaysBitset) = performBFS(history, current.getTaxaBits(), events);

		while (numberOfWaysBitset != 0) {
			int reachable = 63 - __builtin_clzll(numberOfWaysBitset);
			numberOfWaysBitset ^= (1LL << reachable);

			int numberOfWays = numberOfWaysToReach[reachable];

			int startingCount = __builtin_popcount(result.getTaxaBits()) - __builtin_popcount(history);
			int finalCount = __builtin_popcount(result.getTaxaBits()) - __builtin_popcount(reachable);

			double numberOfOptions = getNumberOfOptions(startingCount, finalCount);

			double weight = (double) numberOfWays / numberOfOptions;

			double total = current.getHistory(history) * weight * derivatePuv(startingCount, finalCount, length);

			if (total != 0) {
				result.addToHistory(reachable, total);
			}
		}

	}

	return result;
}

/**
 * Update a list of densemaps.
 */
inline std::vector<densemap> update(const std::vector<densemap>& current, const std::vector<int>& events, double length) {
	std::vector<densemap> result;
	result.reserve(current.size());

	for (auto&& one : current) {
		result.push_back(update(one, events, length));
	}

	return result;
}

/**
 * Update a list of derivates for densemaps.
 */
inline std::vector<densemap> derivativeUpdate(const std::vector<densemap>& current, const std::vector<int>& events, double length) {
	std::vector<densemap> result;
	result.reserve(current.size());

	for (auto&& one : current) {
		result.push_back(derivativeUpdate(one, events, length));
	}

	return result;
}

/**
 * Create every possible subset of the given bitset.
 */
inline std::vector<uint16_t> createSubsets(uint16_t bitset) {
	uint32_t temp = bitset;
	std::vector<uint16_t> result;

	int count = __builtin_popcount(temp);

	result.resize(1 << count);
	for (int i=0;i< (1 << count);i++) {
		result[i] = 0;
	}

	for (int bitNum = 0; bitNum < count; bitNum++) {
		int nextBit = 31 - __builtin_clz(temp);
		temp ^= (1 << nextBit);
		for (int i=0;i< (1 << count);i++) {
			if ((i & (1 << bitNum)) != 0) {
				result[i] |= (1 << nextBit);
			}
		}
	}

	return result;

}

/**
 * Add a result from a split operation.
 */
inline void addResult(const densemap& current, std::vector<densemap>& results, int nodeIndex, uint16_t taxaBits, uint16_t historyBits, int64_t choice, double probability) {
	std::vector<int64_t> choices = current.choices;

	choices[nodeIndex] = choice;

	densemap result;
	result.init(taxaBits, choices);
	result.setHistory(historyBits, probability);

	results.push_back(result);
}

/**
 * Split a densmap at a network node.
 */
inline std::pair<std::vector<densemap>, std::vector<densemap>> split(const std::vector<densemap>& current, int nodeIndex, const std::vector<int>& events, double leftProbability) {
	std::vector<densemap> leftResults;
	std::vector<densemap> rightResults;
	for (auto&& map: current) {
		uint64_t bitset = map.getHistoryBitset();
		while (bitset != 0) {
			int history = 63 - __builtin_clzll(bitset);
			bitset ^= (1LL << history);

			uint16_t fullHistory = map.getTaxaBits() | history;
			for (unsigned int i = 0; i < events.size(); i++) {
				auto& event = events[i];

				// If we had experienced that event
				if (((1 << i) & history) != 0) {
					// We have to remove the corresponding taxabits
					fullHistory &= ~event;
				}
			}

			std::vector<uint16_t> subsets = createSubsets(fullHistory);

			for (unsigned int j = 0; j < subsets.size(); j++) {
				uint16_t finalSubset = subsets[j];

				int numLeft = __builtin_popcount(finalSubset);

				for (unsigned int b = 0; b < events.size(); b++) {
					for (unsigned int i = 0; i < events.size(); i++) {
						auto& event = events[i];

						// If we had experienced that event
						if (((1 << i) & finalSubset) != 0) {

							// We have to add the corresponding taxabits
							finalSubset |= event;
						}
					}
				}
				uint16_t leftSubsetId = j;
				uint16_t rightSubsetId = j ^ (subsets.size() - 1);

				int64_t leftChoiceId =  map.getTaxaBits() | history | (leftSubsetId << 16) | ((long long)fullHistory << 32);
				int64_t rightChoiceId =  map.getTaxaBits() | history | (rightSubsetId << 16) | ((long long)fullHistory << 32);

				uint16_t taxaBits = finalSubset    & 0b1111111111000000;
				uint16_t historyBits = finalSubset & 0b0000000000111111;
				addResult(map, leftResults, nodeIndex, taxaBits, historyBits, leftChoiceId, std::sqrt(map.getHistory(history)) * std::pow(leftProbability, numLeft));
				addResult(map, rightResults, nodeIndex, taxaBits, historyBits, rightChoiceId, std::sqrt(map.getHistory(history)) * std::pow(1-leftProbability, numLeft));
			}
		}
	}

	return { leftResults, rightResults };
}

/**
 * Split the derivatives of a densmap at a network node.
 */
inline std::pair<std::vector<densemap>, std::vector<densemap>> splitDerivatives(const std::vector<densemap>& currentDerivatives, const std::vector<densemap>& current, int nodeIndex, const std::vector<int>& events, double leftProbability) {
	std::vector<densemap> leftResults;
	std::vector<densemap> rightResults;
	for (unsigned int i = 0; i < current.size();i ++) {
		auto&& map = current[i];
		auto&& derivativeMap = currentDerivatives[i];

		uint64_t bitset = map.getHistoryBitset();
		while (bitset != 0) {
			int history = 63 - __builtin_clzll(bitset);
			bitset ^= (1LL << history);

			uint16_t fullHistory = map.getTaxaBits() | history;
			for (unsigned int i = 0; i < events.size(); i++) {
				auto& event = events[i];

				// If we had experienced that event
				if (((1 << i) & history) != 0) {
					// We have to remove the corresponding taxabits
					fullHistory &= ~event;
				}
			}

			std::vector<uint16_t> subsets = createSubsets(fullHistory);

			for (unsigned int j = 0; j < subsets.size(); j++) {
				uint16_t finalSubset = subsets[j];

				int numLeft = __builtin_popcount(finalSubset);

				for (unsigned int b = 0; b < events.size(); b++) {
					for (unsigned int i = 0; i < events.size(); i++) {
						auto& event = events[i];

						// If we had experienced that event
						if (((1 << i) & finalSubset) != 0) {

							// We have to add the corresponding taxabits
							finalSubset |= event;
						}
					}
				}
				uint16_t leftSubsetId = j;
				uint16_t rightSubsetId = j ^ (subsets.size() - 1);

				int64_t leftChoiceId =  map.getTaxaBits() | history | (leftSubsetId << 16) | ((long long)fullHistory << 32);
				int64_t rightChoiceId =  map.getTaxaBits() | history | (rightSubsetId << 16) | ((long long)fullHistory << 32);

				uint16_t taxaBits = finalSubset    & 0b1111111111000000;
				uint16_t historyBits = finalSubset & 0b0000000000111111;
				addResult(map, leftResults, nodeIndex, taxaBits, historyBits, leftChoiceId, derivativeMap.getHistory(history) *  1/(2 * std::sqrt(map.getHistory(history))) * std::pow(leftProbability, numLeft));
				addResult(map, rightResults, nodeIndex, taxaBits, historyBits, rightChoiceId, derivativeMap.getHistory(history) * 1/(2 * std::sqrt(map.getHistory(history))) * std::pow(1-leftProbability, numLeft));
			}
		}
	}

	return { leftResults, rightResults };
}

/**
 * Split a densemap where the derivative is taken with respect to the left probability.
 */
inline std::pair<std::vector<densemap>, std::vector<densemap>> splitDerivativeHere(const std::vector<densemap>& current, int nodeIndex, const std::vector<int>& events, double leftProbability) {
	std::vector<densemap> leftResults;
	std::vector<densemap> rightResults;
	for (auto&& map: current) {
		uint64_t bitset = map.getHistoryBitset();
		while (bitset != 0) {
			int history = 63 - __builtin_clzll(bitset);
			bitset ^= (1LL << history);

			uint16_t fullHistory = map.getTaxaBits() | history;
			for (unsigned int i = 0; i < events.size(); i++) {
				auto& event = events[i];

				// If we had experienced that event
				if (((1 << i) & history) != 0) {
					// We have to remove the corresponding taxabits
					fullHistory &= ~event;
				}
			}

			std::vector<uint16_t> subsets = createSubsets(fullHistory);

			for (unsigned int j = 0; j < subsets.size(); j++) {
				uint16_t finalSubset = subsets[j];

				int numLeft = __builtin_popcount(finalSubset);

				for (unsigned int b = 0; b < events.size(); b++) {
					for (unsigned int i = 0; i < events.size(); i++) {
						auto& event = events[i];

						// If we had experienced that event
						if (((1 << i) & finalSubset) != 0) {

							// We have to add the corresponding taxabits
							finalSubset |= event;
						}
					}
				}
				uint16_t leftSubsetId = j;
				uint16_t rightSubsetId = j ^ (subsets.size() - 1);

				int64_t leftChoiceId =  map.getTaxaBits() | history | (leftSubsetId << 16) | ((long long)fullHistory << 32);
				int64_t rightChoiceId =  map.getTaxaBits() | history | (rightSubsetId << 16) | ((long long)fullHistory << 32);

				uint16_t taxaBits = finalSubset    & 0b1111111111000000;
				uint16_t historyBits = finalSubset & 0b0000000000111111;
				addResult(map, leftResults, nodeIndex, taxaBits, historyBits, leftChoiceId, std::sqrt(map.getHistory(history)) * std::pow(leftProbability, numLeft - 1) * numLeft);
				addResult(map, rightResults, nodeIndex, taxaBits, historyBits, rightChoiceId, std::sqrt(map.getHistory(history)) * std::pow(1-leftProbability, numLeft -1) * -numLeft);
			}
		}
	}

	return { leftResults, rightResults };
}