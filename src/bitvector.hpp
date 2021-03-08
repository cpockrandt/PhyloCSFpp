#pragma once

// for memberset queries in nodes of the newick tree.
// only supports up to 128 species
struct bit_vector {
	uint64_t lo;
	uint64 hi;

	void set_species()
	{
		assert( <= 128);
	}

	bool check_overlap(const bit_vector & other) const noexcept
	{
		return ((other.lo & lo) | (other.hi & hi)) != 0;
	}
}
