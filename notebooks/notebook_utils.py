#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Visualization and plotting utilities for molecular data analysis in Jupyter notebooks.

This module provides reusable functions for exploratory data analysis
of compound and bioactivity data. Functions are designed to work with
traditional scientific modules and are intended to be imported into
Jupyter notebooks.

Available functions:
    - draw_molecule: Draw RDKit molecule object

Usage:
    from notebook_utils import <function_name>
"""

####################################################################################################
## Imports
####################################################################################################

from typing import Any, Dict, List, Tuple

import rdkit
from rdkit.Chem import Draw


####################################################################################################
## Functions
####################################################################################################
def draw_molecule(
    mol: rdkit.Chem.rdchem.Mol,
    legend: str = "",
    width: int = -1,
    height: int = -1,
    highlight_atoms: List[int] | None = None,
    highlight_atom_color: Tuple[float, float, float] | List[Tuple[float, float, float]] = (
        1.0,
        1.0,
        0.7,
    ),
    highlight_bonds: List[int] | None = None,
    highlight_bond_color: Tuple[float, float, float] | List[Tuple[float, float, float]] = (
        0.4,
        0.7,
        1.0,
    ),
    acs_style: bool = False,
) -> Any:
    """Draw molecule and return as SVG text.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object to visualize
    legend : str, optional
        Text label for the molecule, default: ""
    width : int, optional
        Width of the image in pixels, default: -1
    height : int, optional
        Height of the image in pixels, default: -1
    highlight_atoms : List[int] | None , optional
        List of atom indices to highlight, default: None
    highlight_atom_color : Tuple[float, float, float] | List[Tuple[float, float, float]], optional
        Default RGB color for highlighting atoms, default: (1, 1, 0.7)
    highlight_bonds : List[int] | None , optional
        List of bond indices to highlight, default: None
    highlight_bond_color : Tuple[float, float, float] | List[Tuple[float, float, float]], optional
        Default RGB color for highlighting bonds, default: (0.4, 0.7, 1.0)
    acs_style : bool
        Draw molecule with ACS 1996 style, default: False

    Returns
    -------
    str
        SVG image as drawing text
    """
    # Compute 2D coordinates
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    rdkit.Chem.rdDepictor.StraightenDepiction(mol)

    # Create drawer with specified dimensions
    d2d = Draw.MolDraw2DSVG(width, height)

    # Set up atom highlighting
    highlight_atom_map: Dict[int, Tuple[float, float, float]] = {}
    if highlight_atoms is not None:
        if isinstance(highlight_atom_color, tuple):
            highlight_atom_map = dict.fromkeys(highlight_atoms, highlight_atom_color)
        elif isinstance(highlight_atom_color, list):
            highlight_atom_map = dict(zip(highlight_atoms, highlight_atom_color, strict=True))
        else:
            raise ValueError(
                "highlight_atoms must be None or List[int] and highlight_atom_color must"
                + "be Tuple[float, float, float] or List[Tuple[float, float, float]]"
            )

    # Set up bond highlighting
    highlight_bond_map: Dict[int, Tuple[float, float, float]] = {}
    if highlight_bonds is not None:
        if isinstance(highlight_bond_color, tuple):
            highlight_bond_map = dict.fromkeys(highlight_bonds, highlight_bond_color)
        elif isinstance(highlight_bond_color, list):
            highlight_bond_map = dict(zip(highlight_bonds, highlight_bond_color, strict=True))
        else:
            raise ValueError(
                "highlight_bonds must be None or List[int] and highlight_bond_color must"
                + "be Tuple[float, float, float] or List[Tuple[float, float, float]]"
            )

    # Draw the molecule
    draw_options = d2d.drawOptions()  # type: ignore[call-arg]
    if acs_style is not False:
        Draw.SetACS1996Mode(  # type: ignore[attr-defined]
            draw_options,  # type: ignore[call-arg]
            Draw.MeanBondLength(mol),  # type: ignore[attr-defined]
        )
    draw_options.useDefaultAtomPalette()

    d2d.DrawMolecule(
        mol,
        legend=legend,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_map,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_map,
    )
    d2d.FinishDrawing()  # type: ignore[call-arg]

    # Return drawing text
    return d2d.GetDrawingText()  # type: ignore[call-arg]
