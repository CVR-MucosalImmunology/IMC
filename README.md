# Pre-R IMC Pipeline Documentation
## 1. Setting up a `conda` environment

Anaconda is a program used to install packages needed for many steps of the pipeline to run. Follow the steps below to set up Anaconda and a `conda` environment:

**Step 1:** Install [**Anaconda** ](https://www.anaconda.com/download) <br>
**Step 2:** Once Anaconda is installed, navigate to the relevant command line interface:
- **Windows**: 
    - Search for **'Anaconda Prompt'** in the taskbar search
    - Select **Anaconda Prompt**
- **macOS**:
    - Use **`Cmd + Space`** to open Spotlight Search
    - Type **'Terminal'** and press return to open 
    
**Step 3:** Enter the following commands (enter one, press **`Enter`**, then repeat with the next command):
- `git clone --recursive https://github.com/BodenmillerGroup/ImcSegmentationPipeline.git`
- `cd ImcSegmentationPipeline`
- `conda env create -f environment.yml`
- `conda activate imcsegpipe`
- `jupyter lab`

This will automatically open a jupyter instance at `http://localhost:8888/lab` in your browser. From there, open the `1 IMCPreprocessing.ipynb` file and follow the instructions there.

## 2. Installing and opening CellPose

Open **Anaconda Prompt** and enter the following command to install CellPose:
- `conda create -n cellpose pytorch=1.8.2 cudatoolkit=10.2 -c pytorch-lts`

To **open** CellPose (both now and in the future), run both of the following commands:
- `conda activate cellpose`
- `python -m pip install cellpose[gui]` 

## 3. Using the CellPose GUI

**Note:** The steps below were written based on the **CellPose 3** GUI - newer versions may differ slightly

1. Drag an image from the folder into the GUI 
2. Apply the settings below:
- 
3. Optionally set the second channel if you are segmenting cyto and have an available nucleus channel.
4. Click the calibrate button to estimate the size of the objects in the image. Alternatively (RECOMMENDED) you can set the cell diameter by hand and press ENTER. You will see the size you set as a red disk at the bottom left of the image.
5. Click the run segmentation button. If MASKS ON is checked, you should see masks drawn on the image.
Now you can click the LEFT/RIGHT arrow keys to move through the folder and segment another image.