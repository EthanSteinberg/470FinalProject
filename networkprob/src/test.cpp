#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.h"

#include "example.h"

TEST_CASE( "Simple 7 taxa tree case", "[tree]" ) {

	std::vector<TreeNode> genes;
	std::vector<NetNode> species;

	auto prob = calcProbability(createSpecies(species), createGene(genes));

	REQUIRE( prob == Approx(0.000154474) );
}

TEST_CASE( "Simple 7 taxa tree case with trivial introgression", "[treetrivial]" ) {
	std::vector<TreeNode> genes;
	std::vector<NetNode> species;

	auto prob = calcProbability(createSpeciesWithTrivialIntro(species), createGene(genes));

	REQUIRE( prob == Approx(0.000154474) );
}

TEST_CASE( "Simple 7 taxa tree case with fancy introgression", "[introfancy]" ) {
    std::vector<TreeNode> genes;
    std::vector<NetNode> species;

    auto prob = calcProbability(createSpeciesWithIntro(species), createGene(genes));

	REQUIRE( prob == Approx(4.25917e-05) );
}

TEST_CASE( "Simple 7 taxa tree case with fancy introgression by calcDerivatives", "[introfancywithderivative]" ) {
    std::vector<double> derivatives;

    std::vector<TreeNode> genes;
    std::vector<NetNode> species;

    auto prob = calcProbability(createSpeciesWithIntro(species), createGene(genes), &derivatives);

    REQUIRE( prob == Approx(4.25917e-05) );

    REQUIRE( prob == Approx(4.25917e-05) );
}


TEST_CASE( "Make subsets test", "[subset]" ) {

	uint16_t tester = 0b100101;

	std::vector<uint16_t> resultingSets = createSubsets(tester);

	REQUIRE(resultingSets.size() == 8);
	REQUIRE(resultingSets[0] == 0b000000);
	REQUIRE(resultingSets[1] == 0b100000);
	REQUIRE(resultingSets[2] == 0b000100);
	REQUIRE(resultingSets[3] == 0b100100);
	REQUIRE(resultingSets[4] == 0b000001);
	REQUIRE(resultingSets[5] == 0b100001);
	REQUIRE(resultingSets[6] == 0b000101);
	REQUIRE(resultingSets[7] == 0b100101);
}

TEST_CASE( "Test the derivative of puv works properly", "[puvderivative]") {
	double dx = 0.00001;
	double x = 0.1;
	double manual = (puv(3, 2, x + dx) - puv(3, 2, x - dx))/dx * 1.0/2.0;
	double derivative = derivatePuv(3, 2, x);
	REQUIRE(manual == Approx {derivative});
}

TEST_CASE( "Test that split works properly", "[split]" ) {
	densemap source;
	source.init(0b11000000, {-1});
	source.setHistory(0, 1.0);

	std::vector<int> events = { 0b11000000 };

	std::vector<densemap> postUpdate;
	postUpdate.push_back(update(source, events, 1.0));

	REQUIRE(postUpdate[0].getHistory(0) == Approx(0.367879));
	REQUIRE(postUpdate[0].getHistory(1) == Approx(0.632121));

	auto afterSplit = split(postUpdate, 0, events, 0.25);

	REQUIRE(afterSplit.first.size() == 6);
	REQUIRE(afterSplit.second.size() == 6);

	REQUIRE(afterSplit.first[0].getTaxaBits() == 0b00 << 6);
    REQUIRE(afterSplit.first[1].getTaxaBits() == 0b11 << 6);
    REQUIRE(afterSplit.first[2].getTaxaBits() == 0b00 << 6);
	REQUIRE(afterSplit.first[3].getTaxaBits() == 0b10 << 6);
	REQUIRE(afterSplit.first[4].getTaxaBits() == 0b01 << 6);
	REQUIRE(afterSplit.first[5].getTaxaBits() == 0b11 << 6);

    REQUIRE(afterSplit.second[0].getTaxaBits() == 0b00 << 6);
    REQUIRE(afterSplit.second[1].getTaxaBits() == 0b11 << 6);
    REQUIRE(afterSplit.second[2].getTaxaBits() == 0b00 << 6);
    REQUIRE(afterSplit.second[3].getTaxaBits() == 0b10 << 6);
    REQUIRE(afterSplit.second[4].getTaxaBits() == 0b01 << 6);
    REQUIRE(afterSplit.second[5].getTaxaBits() == 0b11 << 6);

    REQUIRE(afterSplit.first[0].getHistory(0) == Approx {0.7950600976});
	REQUIRE(afterSplit.first[0].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.first[1].getHistory(0) == Approx {0});
    REQUIRE(afterSplit.first[1].getHistory(1) == Approx {0.1987650244});

    REQUIRE(afterSplit.first[2].getHistory(0) == Approx {0.6065306597});
    REQUIRE(afterSplit.first[2].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.first[3].getHistory(0) == Approx {0.1516326649});
    REQUIRE(afterSplit.first[3].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.first[4].getHistory(0) == Approx {0.1516326649});
    REQUIRE(afterSplit.first[4].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.first[5].getHistory(0) == Approx {0.0379081662});
    REQUIRE(afterSplit.first[5].getHistory(1) == Approx {0});


    REQUIRE(afterSplit.second[0].getHistory(0) == Approx {0.7950600976});
    REQUIRE(afterSplit.second[0].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.second[1].getHistory(0) == Approx {0});
    REQUIRE(afterSplit.second[1].getHistory(1) == Approx {0.5962950732});

    REQUIRE(afterSplit.second[2].getHistory(0) == Approx {0.6065306597});
    REQUIRE(afterSplit.second[2].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.second[3].getHistory(0) == Approx {0.4548979948});
    REQUIRE(afterSplit.second[3].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.second[4].getHistory(0) == Approx {0.4548979948});
    REQUIRE(afterSplit.second[4].getHistory(1) == Approx {0});

    REQUIRE(afterSplit.second[5].getHistory(0) == Approx {0.3411734961});
    REQUIRE(afterSplit.second[5].getHistory(1) == Approx {0});
}