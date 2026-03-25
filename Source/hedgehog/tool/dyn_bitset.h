#ifndef DYN_BITSET_H
#define DYN_BITSET_H

#include <cstddef>
#include <cstdint>
#include <vector>

/// Dynamically-sized bitset for tracking mesh dependencies.
///
/// Each bit represents a mesh index.  Sized to ceil(N/64) 64-bit words.
/// All operations are O(1) per bit or O(words) for bulk ops.
class DynBitset {
public:
    DynBitset() = default;

    explicit DynBitset(size_t nbits)
        : nbits_(nbits),
          words_((nbits + 63) / 64, 0) {}

    /// Set bit at position `pos`.  Idempotent.
    void set(size_t pos) {
        words_[pos / 64] |= (uint64_t{1} << (pos % 64));
    }

    /// Clear bit at position `pos`.
    void clear(size_t pos) {
        words_[pos / 64] &= ~(uint64_t{1} << (pos % 64));
    }

    /// Test whether bit `pos` is set.
    [[nodiscard]] bool contains(size_t pos) const {
        return (words_[pos / 64] & (uint64_t{1} << (pos % 64))) != 0;
    }

    /// Test whether `other` is a subset of this bitset (all bits in `other`
    /// are also set here).  I.e., `(this & other) == other`.
    [[nodiscard]] bool containsAll(const DynBitset &other) const {
        for (size_t i = 0; i < other.words_.size(); ++i) {
            if ((words_[i] & other.words_[i]) != other.words_[i])
                return false;
        }
        return true;
    }

    /// Number of set bits.
    [[nodiscard]] size_t count() const {
        size_t n = 0;
        for (auto w : words_)
            n += static_cast<size_t>(__builtin_popcountll(w));
        return n;
    }

    /// Total capacity (number of bits).
    [[nodiscard]] size_t size() const { return nbits_; }

    /// Clear all bits.
    void reset() {
        for (auto &w : words_) w = 0;
    }

private:
    size_t nbits_ = 0;
    std::vector<uint64_t> words_;
};

#endif // DYN_BITSET_H
