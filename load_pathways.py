"""
load_pathways.py
================

Load heterologous production pathways into a COBRApy model (e.g. iML1515)
from three CSV files.

Expected schemas
----------------
metabolites.csv : id, name, formula, charge, kegg_id, bigg_id, notes
reactions.csv   : id, name, rxn_str, kegg_id, bigg_id, ec_number, notes
products.csv    : id, name, rxns, product_id, product_tags, precursor_ids, ...

Conventions this loader assumes
-------------------------------
* Metabolite ids in metabolites.csv are given WITHOUT a compartment suffix
  ("nadh"), while the reaction strings reference them WITH a cytosolic suffix
  ("nadh_c"). Every heterologous species is placed in the cytosol ("_c").
  If your base model uses a different id for a native metabolite, the loader
  matches on "<id>_c" and reuses the model's copy rather than duplicating it.
* "=>" is treated as irreversible (forward), "<=>" as reversible, exactly as
  cobra's build_reaction_from_string parses them.

Typical use in a notebook
--------------------------
    import cobra
    from load_pathways import build_extended_model

    model, products = build_extended_model(
        base_model="iML1515",              # BiGG id, or a path to an SBML file
        metabolites_csv="metabolites.csv",
        reactions_csv="reactions.csv",
        products_csv="products.csv",
    )
"""

import ast
import cobra
import pandas as pd
from cobra import Metabolite, Reaction

DEFAULT_COMPARTMENT = "c"


# --------------------------------------------------------------------------- #
# small parsing helpers
# --------------------------------------------------------------------------- #
def _clean(value):
    """Return a stripped string, or None for NaN/empty."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip()
    return s or None


def _met_id(bare_id, compartment=DEFAULT_COMPARTMENT):
    """metabolites.csv id -> in-model id, e.g. 'nadh' -> 'nadh_c'."""
    return f"{bare_id}_{compartment}"


def _parse_list(cell):
    """
    Parse a cell that may be a Python list literal ("['a', 'b']"), a bare token
    ("accoa"), or empty. Always returns a list of strings.
    """
    s = _clean(cell)
    if s is None:
        return []
    try:
        val = ast.literal_eval(s)
    except (ValueError, SyntaxError):
        return [s]
    if isinstance(val, (list, tuple)):
        return [str(x).strip() for x in val]
    return [str(val).strip()]


# --------------------------------------------------------------------------- #
# metabolites
# --------------------------------------------------------------------------- #
def load_metabolites(model, met_df, compartment=DEFAULT_COMPARTMENT,
                     overwrite_existing=False, verbose=True):
    """
    Add each metabolite in `met_df` to `model` in the given compartment.

    A metabolite already present in the model (matched on '<id>_<compartment>')
    is left untouched unless `overwrite_existing=True`, in which case its
    formula/charge are updated from the CSV. Returns (added_ids, skipped_ids).
    """
    added, skipped = [], []
    for _, row in met_df.iterrows():
        bare = _clean(row["id"])
        if bare is None:
            continue
        mid = _met_id(bare, compartment)
        formula = _clean(row.get("formula"))
        charge = _clean(row.get("charge"))
        charge = int(float(charge)) if charge is not None else None

        if mid in model.metabolites:
            if overwrite_existing:
                m = model.metabolites.get_by_id(mid)
                if formula is not None:
                    m.formula = formula
                if charge is not None:
                    m.charge = charge
            skipped.append(mid)
            continue

        met = Metabolite(
            mid,
            name=_clean(row.get("name")) or "",
            formula=formula,
            charge=charge,
            compartment=compartment,
        )
        kegg = _clean(row.get("kegg_id"))
        if kegg:
            met.annotation["kegg.compound"] = kegg
        model.add_metabolites([met])
        added.append(mid)

    if verbose:
        print(f"[metabolites] added {len(added)}, reused {len(skipped)} existing")
    return added, skipped


# --------------------------------------------------------------------------- #
# reactions
# --------------------------------------------------------------------------- #
def load_reactions(model, rxn_df, default_bounds=(0.0, 1000.0),
                   skip_existing=True, verbose=True):
    """
    Add each reaction in `rxn_df` to `model`, building stoichiometry from the
    `rxn_str` column. Direction ("=>" vs "<=>") is honoured by cobra's parser.

    Returns (added_ids, skipped_ids). Any metabolite referenced in a reaction
    string but not already in the model is auto-created by cobra WITHOUT a
    formula; those are reported so you can catch typos early.
    """
    added, skipped = [], []
    mets_before = {m.id for m in model.metabolites}

    for _, row in rxn_df.iterrows():
        rid = _clean(row["id"])
        rxn_str = _clean(row["rxn_str"])
        if rid is None or rxn_str is None:
            continue
        if rid in model.reactions:
            if skip_existing:
                skipped.append(rid)
                continue

        rxn = Reaction(rid)
        rxn.name = _clean(row.get("name")) or ""
        rxn.lower_bound, rxn.upper_bound = default_bounds  # overwritten by parser
        model.add_reactions([rxn])
        # build_reaction_from_string must run after the reaction is in a model
        rxn.build_reaction_from_string(rxn_str)

        ec = _clean(row.get("ec_number"))
        if ec:
            rxn.annotation["ec-code"] = ec
        added.append(rid)

    # anything cobra had to invent while parsing (usually a typo)
    orphans = [m.id for m in model.metabolites
               if m.id not in mets_before and m.formula is None]
    if verbose:
        print(f"[reactions]   added {len(added)}, skipped {len(skipped)} existing")
        if orphans:
            print(f"[reactions]   WARNING auto-created metabolites with no formula "
                  f"(check for typos): {sorted(orphans)}")
    return added, skipped


# --------------------------------------------------------------------------- #
# products
# --------------------------------------------------------------------------- #
def load_products(products_df, compartment=DEFAULT_COMPARTMENT):
    """
    Parse products.csv into a dict keyed by module id (the `id` column, e.g.
    'etoh_pdc'). Each value carries the module's reaction list, the target
    metabolite id in the model, tags, precursors and formula.
    """
    products = {}
    for _, row in products_df.iterrows():
        pid = _clean(row["id"])
        if pid is None:
            continue
        target = _clean(row.get("product_id"))
        products[pid] = {
            "name": _clean(row.get("name")),
            "reactions": _parse_list(row.get("rxns")),
            "product_id": target,
            "product_met": _met_id(target, compartment) if target else None,
            "tags": _parse_list(row.get("product_tags")),
            "precursors": _parse_list(row.get("precursor_ids")),
            "formula": _clean(row.get("formula")),
        }
    return products


def add_product_boundaries(model, products, boundary_type="demand", verbose=True):
    """
    Give each product metabolite a way to leave the system so it can carry flux.

    NOTE ON MODELLING CHOICE: with no curated transport/secretion reactions in
    the CSVs, this adds a cytosolic drain per product. type="demand" -> DM_<met>
    (irreversible, lb=0), type="sink" -> reversible, type="exchange" -> EX_<met>.
    A demand reaction on the cytosolic pool is the usual simplification for
    yield/production-envelope work, but confirm it matches the convention Galib
    /Dr. Trinh want before you report numbers -- it changes the accounting.

    Returns {module_id: boundary_reaction_id}. Reuses an existing boundary if the
    product metabolite already has one.
    """
    created = {}
    for pid, info in products.items():
        mid = info["product_met"]
        if mid is None or mid not in model.metabolites:
            if verbose:
                print(f"[boundaries]  skip {pid}: metabolite {mid} not in model")
            continue
        met = model.metabolites.get_by_id(mid)
        existing = [r for r in met.reactions if r.boundary]
        if existing:
            created[pid] = existing[0].id
            continue
        rxn = model.add_boundary(met, type=boundary_type)
        created[pid] = rxn.id
    if verbose:
        print(f"[boundaries]  ensured {len(created)} product drains "
              f"(type='{boundary_type}')")
    return created


# --------------------------------------------------------------------------- #
# validation helpers
# --------------------------------------------------------------------------- #
def check_pathway_balance(model, rxn_df, verbose=True):
    """
    Run cobra's check_mass_balance on every reaction from rxn_df that made it
    into the model. Returns {rxn_id: {element/charge: imbalance}} for anything
    that does not close. An empty dict means every heterologous reaction is
    mass- and charge-balanced (given the formulas/charges in the model).
    """
    imbalanced = {}
    for rid in rxn_df["id"].map(_clean).dropna():
        if rid not in model.reactions:
            continue
        try:
            imb = model.reactions.get_by_id(rid).check_mass_balance()
        except Exception as exc:  # e.g. a participating metabolite has no formula
            imb = {"error": str(exc)}
        if imb:
            imbalanced[rid] = imb
    if verbose:
        if imbalanced:
            print(f"[balance]     {len(imbalanced)} reaction(s) NOT balanced:")
            for rid, imb in imbalanced.items():
                print(f"                {rid}: {imb}")
        else:
            print("[balance]     all heterologous reactions balanced")
    return imbalanced


def validate_product_reactions(products, rxn_df, model=None):
    """
    Confirm every reaction id referenced by a product module exists either in
    reactions.csv or (if `model` given) already in the base model. Returns
    {module_id: [missing_reaction_ids]} for modules with gaps.
    """
    known = set(rxn_df["id"].map(_clean).dropna())
    if model is not None:
        known |= {r.id for r in model.reactions}
    missing = {}
    for pid, info in products.items():
        gaps = [r for r in info["reactions"] if r not in known]
        if gaps:
            missing[pid] = gaps
    return missing


def get_module_reactions(products, module_id):
    """Reaction ids that make up one product module."""
    return list(products[module_id]["reactions"])


# --------------------------------------------------------------------------- #
# orchestrator
# --------------------------------------------------------------------------- #
def build_extended_model(base_model, metabolites_csv, reactions_csv,
                         products_csv=None, compartment=DEFAULT_COMPARTMENT,
                         add_boundaries=True, boundary_type="demand",
                         check_balance=True, verbose=True):
    """
    Load a base model, add all heterologous metabolites + reactions, optionally
    add product drains, and validate. Returns (model, products_dict).

    `base_model` may be:
        * a cobra.Model            -> copied, original untouched
        * a BiGG model id string   -> fetched via cobra.io.load_model (needs net)
        * a path to an SBML/JSON   -> read from disk
    """
    if isinstance(base_model, cobra.Model):
        model = base_model.copy()
    elif isinstance(base_model, str):
        try:
            model = cobra.io.load_model(base_model)
        except Exception:
            model = cobra.io.read_sbml_model(base_model)
    else:
        raise TypeError("base_model must be a cobra.Model, BiGG id, or file path")

    met_df = pd.read_csv(metabolites_csv)
    rxn_df = pd.read_csv(reactions_csv)

    load_metabolites(model, met_df, compartment, verbose=verbose)
    load_reactions(model, rxn_df, verbose=verbose)

    products = {}
    if products_csv is not None:
        products = load_products(pd.read_csv(products_csv), compartment)
        gaps = validate_product_reactions(products, rxn_df, model)
        if gaps and verbose:
            print(f"[products]    modules referencing unknown reactions: {gaps}")
        if add_boundaries:
            add_product_boundaries(model, products, boundary_type, verbose=verbose)

    if check_balance:
        check_pathway_balance(model, rxn_df, verbose=verbose)

    return model, products


if __name__ == "__main__":
    # quick self-test against a fresh empty model (no base model / no network):
    # builds a self-contained model from the CSVs and checks every reaction's
    # mass/charge balance using the formulas provided.
    m = cobra.Model("heterologous_only")
    mets = pd.read_csv("metabolites.csv")
    rxns = pd.read_csv("reactions.csv")
    load_metabolites(m, mets)
    load_reactions(m, rxns)
    prods = load_products(pd.read_csv("products.csv"))
    print("missing product reactions:", validate_product_reactions(prods, rxns))
    check_pathway_balance(m, rxns)
