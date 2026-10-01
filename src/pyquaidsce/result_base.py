"""Shared named-coefficient and elasticity output for demand systems."""

from typing import Dict


class DemandResults:
    def named(self) -> Dict[str, float]:
        return dict(zip(self.names, self.b))

    def get(self, key: str) -> float:
        """Look up a coefficient by ``eq:name`` or bare ``name``."""
        d = self.named()
        if key in d:
            return d[key]
        for k, v in d.items():
            if k.split(":", 1)[-1] == key:
                return v
        raise KeyError(key)

    # ------------------------------------------------------------------ #
    def elasticity_tables(self) -> str:
        n = self.spec.neqn
        nm = self.share_names
        out = ["", "Expenditure (income) elasticities, at means", "-" * 44]
        w = max(len(x) for x in nm) + 2
        for i in range(n):
            out.append(f"  {nm[i]:<{w}} {self.elas.income[i]:>12.6f}")
        for lab, M in (
            ("Uncompensated (Marshallian) price elasticities [row = good, "
             "column = price]", self.elas.uncompensated),
            (getattr(self, "compensated_label", "Compensated (Hicksian) price elasticities")
             + " [row = good, column = price]", self.elas.compensated),
        ):
            out += ["", lab, "-" * min(len(lab), 100)]
            out.append(" " * w + "".join(f"{x[:10]:>12}" for x in nm))
            for i in range(n):
                out.append(
                    f"  {nm[i]:<{w}}" + "".join(f"{M[i, j]:>12.6f}"
                                                for j in range(n))
                )
        return "\n".join(out)
