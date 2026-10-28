# Thesis supplementary code

ImageJ/Fiji macros and the geometric vibratome-section model used for my MSc thesis. This repository contains analysis code only.

The interactive explorer is the same program as [Bram6921/Vibratome_sections](https://github.com/Bram6921/Vibratome_sections). The figure script below is included here because it is the script that generated the thesis figure.

## Layout

| Path | Role |
| --- | --- |
| `imagej/follicle_only_measurement.ijm` | Follicle area, major axis, minor axis, and eccentricity |
| `imagej/Distance_measurements.ijm` | Cortical depth from a straight line, calibrated from acquisition XML |
| `imagej/Cryo_Confocal_Target_Screening.ijm` | Mark targets on a cryo-confocal CZI stack |
| `imagej/Cryo_Confocal_Collaborator_Export.ijm` | Export marked overviews, region crops, and target zooms from a screening result |
| `geometric_model/geometry.py` | Volume engine: knife slab, tissue block, and intersection volumes |
| `geometric_model/app.py` | Browser explorer. Also imported by the figure script for the 3D panels |
| `geometric_model/make_parameter_sweep_geometric_model_final.py` | Parameter sweeps and the thesis figure |
| `requirements.txt` | Python packages for the geometric model |

## ImageJ macros

Install [Fiji](https://fiji.sc/) with the Bio-Formats plugin. The cryo-confocal macros were written for Fiji / ImageJ 1.54. The follicle and depth macros require ImageJ 1.53 and 1.52p.

Run a macro with **Plugins → Macros → Run…** and select the `.ijm` file. Do not press Cancel on a macro prompt; Cancel stops the macro.

The macros ask for an input folder and an output folder when they start. They do not read a fixed path on disk.

### Follicle morphometry

`follicle_only_measurement.ijm`

1. Choose the experiment folder. It should contain one subfolder per z-stack.
2. Each stack is a TIFF series. The macro opens the file ending in `_z00.tif` or `_z00.tiff` and groups the series with Bio-Formats.
3. For automatic pixel size and z-step, put the acquisition file at `MetaData/<stackBase>_Properties.xml` inside that stack folder. If it is missing, the macro asks for the pixel size and z-step.
4. For each follicle: mark the bottom and top planes and whether each pole is visible, draw one ROI on the clearest slice, then classify follicle type and granulosa shape.
5. Measurements from every stack are appended to one semicolon-separated CSV. ROI files are saved beside it.

Recorded quantities are area, major axis, minor axis, and eccentricity from the fitted ellipse. Z-extent is calculated only when both poles are visible.

### Cortical depth

`Distance_measurements.ijm`

1. Choose the folder that contains the TIFF file or TIFF sequence and a `*_Properties.xml` file in that folder or one subfolder down.
2. The macro reads the X voxel size from that XML and sets the image scale.
3. Draw a straight line from the ovarian surface to a follicle. Choose `primary follicle`, `primordial follicle`, or `not distinguishable`.
4. Repeat until finished, then choose the CSV filename.

The CSV columns are measurement number, calibrated distance, unit, pixel size, and follicle type.

### Target screening

`Cryo_Confocal_Target_Screening.ijm`

1. Choose an output folder, then a `.czi` z-stack. Multiple channels are supported. There should be one time point.
2. Scroll in Z, zoom, and draw an oval around each target.
3. OK with an oval selected adds that target. OK without a selection finishes the run.
4. Each target stays visible on every Z slice. The label records the slice where it was drawn, for example `1 (z: 4/9)`.

Saved files, using the CZI basename:

- `<name>.czi`, a copy of the input
- `<name>_targets.zip`, the ROI set
- `<name>_targets.csv`, pixel coordinates and the Z slice of each target
- `<name>_max_projection.png`, a max projection with the marks burned in

### Collaborator export

`Cryo_Confocal_Collaborator_Export.ijm`

Run this after target screening. Put `<name>.czi` and `<name>_targets.zip` in one folder.

The macro writes:

- `<name>_max_confocal.png`, confocal max projection with marks
- `<name>_max_reflection.png`, reflection max projection with marks, after a brightness and contrast adjustment
- `<name>_intermediate_crop_XX_confocal.png` and `_reflection.png`, 2000×2000 px crops of a region you place
- `<name>_intermediate_crops.csv`, crop coordinates
- `<name>_zoomed_targetXX_zZZ_confocal.png`, a 350×350 px confocal crop on the Z slice where that target was marked, with a 50 µm scale bar and no marks

The scale bar uses a pixel size you enter when the macro asks, for example 100 px = 10 µm.

Screening is the step that records targets. The export macro only makes images from those saved ROIs.

## Geometric model

Requires Python 3.10 or newer.

From this repository:

```bash
python -m pip install -r requirements.txt
cd geometric_model
python app.py
```

Then open http://127.0.0.1:8050. Leave the terminal open. Closing it stops the explorer.

If `python` is not found, use `py` or `python3`.

### What the model calculates

The tissue block starts as a 3 × 3 × 3 mm cube. A target layer of 3 × 3 × 0.2 mm sits 0.1 mm below the top face. A section is a 0.1 mm knife slab, infinite in the cutting plane and clipped to the tissue. Volumes are the intersection of that slab with the tissue and with the target. Empty space outside the block is not part of a section.

The top face can be flat or curved. The target can stay a cuboid or follow the top surface.

Section numbers follow collection order. Section 1 is the first 0.1 mm slab whose outer boundary touches the surface-facing support plane. Later sections move into the tissue.

### Thesis figure

From `geometric_model`, after the packages above are installed:

```bash
python make_parameter_sweep_geometric_model_final.py
```

The script imports `geometry.py` and `app.py` in this folder. It writes CSV, JSON, PNG, and PDF files to `geometric_model/output/`. The figure file is `output/parameter_sweep_geometric_model_V6.pdf`. The run computes several orientation and curvature sweeps on a 128×128 grid, so it can take a while.

Those outputs are generated files and are not stored in this repository.
