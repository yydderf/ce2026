from abc import ABC, abstractmethod
from fractions import Fraction

import math
import random
import argparse


# ref: https://en.wikipedia.org/wiki/Curve25519

def modular_inverse(x, P):
    # modular inverse via Fermat's Little Theorem (x^{p-1} == 1 mod P)
    return pow(x, P - 2, P)

# the square root of -1 (used for point decompression)
# I^2 == pow(2, (P-1) // 2, P) == P - 1 by Euler's criterion
# I = pow(2, (P - 1) // 4, P)

# only abstract functions that is shared among all possible schemes (e.g. ed448, ed25519, etc..)
class TableBuilder:
    def __init__(self, window_size: int = 4, scalar_bits: int = 256):
        if not isinstance(window_size, int) or window_size < 1:
            raise ValueError("window_size must be a positive integer")
        self.window_size = window_size
        self.scalar_bits = scalar_bits
        self._table = None

    @property
    def radix(self) -> int:
        return 1 << self.window_size

    @property
    def entries_per_table(self) -> int:
        return self.radix // 2    

    @property
    def num_tables(self) -> int:
        # ceil(scalar_bits / w) digits, plus one for the carry of signed recoding
        return -(-self.scalar_bits // self.window_size) + 1

    # curve specific ---
    @property
    @abstractmethod
    def identity(self):
        ...
 
    @property
    @abstractmethod
    def base_point(self):
        ...
 
    @abstractmethod
    def point_add(self, p, q):
        ...
 
    @abstractmethod
    def point_double(self, p):
        ...
 
    @abstractmethod
    def point_neg(self, p):
        ...
 
    @abstractmethod
    def to_affine(self, p):
        ...

    # io
    @abstractmethod
    def to_c_header(self, *args, **kwargs) -> str:
        ...

    def build_table(self):
        """Compute (and cache) the table. Evens by doubling, odds by adding."""
        if self._table is not None:
            return self
        self._table, cur = [], self.base_point          # cur = radix^i * B
        for _ in range(self.num_tables):
            row = [None, cur]
            for j in range(2, self.entries_per_table + 1):
                row.append(self.point_double(row[j // 2]) if j % 2 == 0
                           else self.point_add(row[j - 1], cur))
            self._table.append(row[1:])
            for _ in range(self.window_size):     # cur *= radix  (w doublings)
                cur = self.point_double(cur)
        return self
 
    def recode_signed(self, k: int):
        """k = sum d_i * radix^i, with d_i in [-(radix/2 - 1), radix/2]."""
        digits, half = [], self.radix // 2
        for _ in range(self.num_tables):
            d = k & (self.radix - 1)
            k >>= self.window_size
            if d > half:
                d -= self.radix
                k += 1
            digits.append(d)
        if k != 0:
            raise ValueError("scalar too large for this table")
        return digits
 
    def fixed_base_mult(self, k: int):
        """Reference k*B using the table. Variable-time: for testing only."""
        table, r = self.build_table(), self.identity
        for i, d in enumerate(self.recode_signed(k)):
            if d > 0:
                r = self.point_add(r, table[i][d - 1])
            elif d < 0:
                r = self.point_add(r, self.point_neg(table[i][-d - 1]))
        return r
 
    def scalar_mult(self, k: int, p):
        """Plain double-and-add, independent of the table (for cross-checking)."""
        r = self.identity
        for bit in bin(k)[2:]:
            r = self.point_double(r)
            if bit == "1":
                r = self.point_add(r, p)
        return r


# inherit from table builder
class Ed25519TableBuilder(TableBuilder):
    """Ed25519: -x^2 + y^2 = 1 + d x^2 y^2 over GF(2^255 - 19), extended coordinates."""
    FIELD_BITS = 255
    P = 2 ** 255 - 19
    # curve coefficient
    # pow(121666, P - 2, P) % P -> modular inverse via Fermat's Little Theorem, 121666^-1
    # related to the Montgomery form's A = 486662, (A+2)/4 = 121666
    D = (-121665 * modular_inverse(121666, P)) % P
    # the large prime factored from #E
    L = 2 ** 252 + 27742317777372353535851937790883648493
    BY = (4 * modular_inverse(5, P)) % P
    BX = 15112221349535400772501151409588531511454012693041857206046113283949847762202

    def __init__(self, window_size=4, scalar_bits=256):
        super().__init__(window_size, scalar_bits)

    # points are (X, Y, Z, T), x = X/Z, y = Y/Z, xy = T/Z
    @property
    def identity(self):
        return (0, 1, 1, 0)
 
    @property
    def base_point(self):
        return (self.BX, self.BY, 1, self.BX * self.BY % self.P)
 
    def point_add(self, p1, p2):
        P = self.P
        X1, Y1, Z1, T1 = p1
        X2, Y2, Z2, T2 = p2
        A = (Y1 - X1) * (Y2 - X2) % P
        B = (Y1 + X1) * (Y2 + X2) % P
        C = 2 * T1 * T2 * self.D % P
        D_ = 2 * Z1 * Z2 % P
        E, F, G, H = (B - A) % P, (D_ - C) % P, (D_ + C) % P, (B + A) % P
        return E * F % P, G * H % P, F * G % P, E * H % P
 
    def point_double(self, p):
        P = self.P
        X1, Y1, Z1, _ = p
        A, B, C = X1 * X1 % P, Y1 * Y1 % P, 2 * Z1 * Z1 % P
        # E, F, G, H are negated w.r.t. the textbook formulas; the signs cancel in every product below
        H = (A + B) % P
        E = (H - (X1 + Y1) ** 2) % P
        G = (A - B) % P
        F = (C + G) % P
        return E * F % P, G * H % P, F * G % P, E * H % P
 
    def point_neg(self, p):
        X, Y, Z, T = p
        return (-X % self.P, Y, Z, -T % self.P)
 
    def to_affine(self, p):
        X, Y, Z, _ = p
        iz = modular_inverse(Z, self.P)
        return X * iz % self.P, Y * iz % self.P

    def limb_layout(self, radix: float):
        e = Fraction(radix).limit_denominator(1000)
        n = Fraction(self.FIELD_BITS) / e
        if e <= 0 or n.denominator != 1:
            raise ValueError(f"limbs of radix 2^{radix} do not tile {self.FIELD_BITS} bits")
        return [math.ceil(e * i) for i in range(int(n) + 1)]
 
    def to_limbs(self, v: int, starts):
        return [(v >> starts[i]) & ((1 << (starts[i + 1] - starts[i])) - 1) for i in range(len(starts) - 1)]
 
    def to_c_header(self, radix: float = 25.5) -> str:
        """Emit the table as (y+x, y-x, 2d*x*y) entries split into 2^radix limbs.
 
        radix sets the limb size: 25.5 -> 10 x int32_t, 51 -> 5 x int64_t.
        """
        if self._table is None:
            raise ValueError("table has not been built yet")
        starts = self.limb_layout(radix)
        n_limbs = len(starts) - 1
        widest = max(b - a for a, b in zip(starts, starts[1:]))
        if widest <= 31:
            limb_t, sfx = "int32_t", ""
        elif widest <= 63:
            limb_t, sfx = "int64_t", "LL"
        else:
            raise ValueError("limbs wider than 63 bits are not supported")
 
        def fe(v):
            return "{" + ", ".join(f"{x}{sfx}" for x in self.to_limbs(v, starts)) + "}"
 
        def entry(p):
            x, y = self.to_affine(p)
            P = self.P
            return "{" + ", ".join(fe(v) for v in ((y + x) % P, (y - x) % P, 2 * self.D * x * y % P)) + "}"
 
        guard = f"ED25519_BASE_TABLE_H"
        out = [
            f"/* Ed25519 fixed-base table: window {self.window_size}, {self.num_tables} tables x {self.entries_per_table} entries.",
            f"   base_table[i][j-1] = affine (x || y) of j * 2^({self.window_size}*i) * B,",
            "   each coordinate 32 bytes little-endian. Generated file. */",
            f"#ifndef {guard}",
            f"#define {guard}",
            "",
            f"#define ED25519_BASE_WINDOW  {self.window_size}",
            f"#define ED25519_BASE_RADIX   {1 << self.window_size}",
            f"#define ED25519_BASE_TABLES  {self.num_tables}",
            f"#define ED25519_BASE_NDIGITS {self.num_tables - 1}",
            f"#define ED25519_BASE_ENTRIES {self.entries_per_table}",
            f"#define ED25519_FE_LIMBS     {n_limbs}",
            "",
            f"static const unsigned char ed25519_base_table[{self.num_tables}][{self.entries_per_table}][64] = {{",
        ]
        for i, row in enumerate(self._table):
            out.append(f"  {{ /* table {i} */")
            for p in row:
                x, y = self.to_affine(p)
                raw = x.to_bytes(32, "little") + y.to_bytes(32, "little")
                out.append("    {" + ",".join(f"0x{b:02x}" for b in raw) + "},")
            out.append("  },")
        out += ["};", "", f"#endif /* {guard} */", ""]
        return "\n".join(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an Ed25519 fixed-base table header.")
    parser.add_argument("--window", type=int, default=4, help="window size w in bits")
    parser.add_argument("--limb-radix", type=float, default=25.5, help="limb radix exponent: 25.5 (32-bit) or 51 (64-bit)")
    parser.add_argument("--out", default=None, help="output header path")
    args = parser.parse_args()
 
    out = args.out or f"ed25519_base_w{args.window}_r{str(args.limb_radix).replace('.', '_')}.h"
    with open(out, "w") as f:
        f.write(Ed25519TableBuilder(window_size=args.window).build_table().to_c_header(radix=args.limb_radix))
    print("done")
