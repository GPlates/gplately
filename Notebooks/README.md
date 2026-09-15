# GPlately Notebooks

**📌 Please remove the outputs of all cells before checking in notebooks.** The images in the
outputs take up a large amount of space and cause GitHub repository bloat. Visit
[this page](https://gplates.github.io/gplately/latest/sphinx/html/examples.html) if you would
like to see the outputs.

## Running the notebooks

You can run these notebooks either in a local micromamba/conda environment or with the GPlately
Docker image.

### Option 1: micromamba environment

Create and activate an environment from [`docker/env.yaml`](../docker/env.yaml) — this is the
same environment used to build the Docker image, and already includes Jupyter, PyGMT and the
other packages the notebooks need — then launch Jupyter from this `Notebooks` folder:

```bash
# from the repository root
micromamba create -f docker/env.yaml
micromamba activate gplately-env-for-docker
cd Notebooks
jupyter notebook
```

### Option 2: GPlately Docker image

Pull the image and run it, mounting a local directory so any files you create are kept on your
host machine:

```bash
docker pull gplates/gplately
docker run --rm -ti -p 8888:8888 -v `pwd`:/workspace/my_stuff gplately
```

Then open the URL Jupyter prints in the terminal (`http://127.0.0.1:8888/...`) in your browser —
the notebooks are already copied into the image, under `/workspace/Notebooks`. See
[`docker/README.md`](../docker/README.md) for building the image yourself and other Docker
options.

📌 Note: [`12-MutschlerWorldPorphyryCopperDepositsRegionalPlots.ipynb`](12-MutschlerWorldPorphyryCopperDepositsRegionalPlots.ipynb)
needs extra data downloaded before it will run — it will not work out of the box.

## Basic Examples

The [`Examples`](Examples) folder contains short, focused notebooks that each demonstrate a
single GPlately task:

- [**01 - Hello World**](Examples/01-HelloWorld.ipynb): The minimal working example — create a paleo-map with GPlately in a few lines.
- [**02 - Plate Model Manager**](Examples/02-PlateModelManager.ipynb): Use `PlateModelManager` to access plate model files.
- [**03 - Plot with Cartopy**](Examples/03-PlotWithCartopy.ipynb): Build a more detailed paleo-map with Cartopy.
- [**04 - Plot with PyGMT**](Examples/04-PlotWithPyGMT.ipynb): Plot maps using GPlately's PyGMT integration.
- [**05 - Reconstruct Files**](Examples/05-ReconstructFiles.ipynb): Reconstruct shapefiles with GPlately.
- [**06 - Load Plate Model from Files**](Examples/06-LoadPlateModelFromFiles.ipynb): Use your own local plate model and reconstruct points.
- [**07 - Save Reconstructed Data**](Examples/07-SaveReconstructedData.ipynb): Save reconstructed data to shapefiles and other formats.
- [**08 - Use Auxiliary Functions**](Examples/08-UseAuxiliaryFunctions.ipynb): Quickly create `PlateReconstruction` and `PlotTopologies` objects with `gplately.auxiliary`.
- [**09 - Icosahedron Mesh**](Examples/09-IcosahedronMesh.ipynb): Generate and plot an icosahedron mesh.
- [**10 - Color Map and Color Palette Table**](Examples/10-ColorMapAndColorPaletteTable.ipynb): Work with Matplotlib colormaps and GMT Color Palette Tables (CPT).

(The `Examples` folder also contains `PNG_reconstruction_copper_deposits.ipynb`, a larger
standalone worked example, which is not listed above.)

## Sample workflows

To see GPlately in action, launch a Jupyter Notebook environment and check out the sample
notebooks listed below.

- [**01 - Getting Started**](01-GettingStarted.ipynb): A brief overview of how to initialise GPlately's main objects.
- [**02 - Plate Reconstructions**](02-PlateReconstructions.ipynb): Setting up a `PlateReconstruction` object, reconstructing geological data through time.
- [**03 - Working with Points**](03-WorkingWithPoints.ipynb): Setting up a `Points` object, reconstructing seed point locations through time. This notebook uses point data from the Paleobiology Database (PBDB).
- [**04 - Velocity Basics**](04-VelocityBasics.ipynb): Calculating plate velocities, plotting velocity vector fields.
- [**05 - Working with Feature Geometries**](05-WorkingWithFeatureGeometries.ipynb): Processing and plotting assorted polyline, polygon and point data from [GPlates 2.3's sample data sets](https://www.earthbyte.org/gplates-2-3-software-and-data-sets/).
- [**06 - Rasters**](06-Rasters.ipynb): Reading, resizing and resampling raster data, and linearly interpolating point data onto raster data.
- [**07 - Plate Tectonic Stats**](07-WorkingWithPlateTectonicStats.ipynb): Calculate and plot subduction zone and ridge data (convergence/spreading velocities, subduction angles, subduction zone and ridge lengths, crustal surface areas produced and subducted, etc.).
- [**08 - Predicting Slab Flux**](08-PredictingSlabFlux.ipynb): Predicting the average slab dip angle of subducting oceanic lithosphere.
- [**09 - Motion Paths and Flowlines**](09-CreatingMotionPathsAndFlowlines.ipynb): Using pyGPlates to create motion paths and flowlines of points on a tectonic plate to illustrate the plate's trajectory through geological time.
- [**10 - Seafloor Grid**](10-SeafloorGrids.ipynb): Defines the parameters needed to set up a `SeafloorGrid` object, and demonstrates how to produce age and spreading rate grids from a set of plate reconstruction model files.
- [**11 - Andes Fluxes**](11-AndesFluxes.ipynb): Demonstrates how the reconstructed subduction history along the Andean margin can be used in plate kinematics analysis and data mining.
- [**12 - Mutschler World Porphyry Copper Deposits Regional Plots**](12-MutschlerWorldPorphyryCopperDepositsRegionalPlots.ipynb): Generates regional plots for Mutschler world porphyry copper deposits.
- [**13 - Reconstructing Zircon Data**](13-ReconstructingZirconData.ipynb): Demonstrates how to reconstruct and plot zircon data on a global map through geological time.
- [**14 - Rule Based GPML Processing Pipeline**](14-RuleBasedGPMLProcessingPipeline.ipynb): Demonstrates how to use the rule-based GPML processing pipeline to filter and transform geological data.
- [**15 - Convert Grid Reference Frame**](15-ConvertGridReferenceFrame.ipynb): Demonstrates how to convert the reference frame of grids.
