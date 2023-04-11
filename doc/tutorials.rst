Tutorials
=========

This page contains two tutorials:

1. :ref:`Read and explore OCT data <tutorial_read_oct>`:
2. :ref:`Feature extraction <tutorial_feature_extraction>`

.. _tutorial_read_oct:

**Read and explore OCT data**
-----------------------------
This tutorial shows how to:

* Load OCT files into MATLAB
* Visualize loaded data (fundus image, OCT slices, segmentation...)
* Generate a visual report

The source code of this tutorial can be found in the ``tutorials`` folder or `here <https://github.com/drombas/retimat/blob/main/tutorials/File_Reading_Visualization.mlx>`_.

0. Understand your data
^^^^^^^^^^^^^^^^^^^^^^^
The first step is to understand your images in terms of:

* *Device and file format*: there are many different devices and file formats and not all of them are equally easy to read and have the same information.
* *Acquisition protocol*: OCT images are composed of one or more 2D images (aka slices or B-scans). It is crucial to know the area covered by these B-scans (macula, optic disk, wide), as well as their geometrical arrangement (horizontal, circular).

Understanding the above will help troubleshoots problems.

RETIMAT can only read vol, e2e, fda, img files. If your files are different, you can try the following Python packages:

* `OCT-Converter <hhttps://github.com/marksgraham/OCT-Converter>`_ 
* `eyepy <https://github.com/MedVisBonn/eyepy>`_ 

1. Load the images
^^^^^^^^^^^^^^^^^^
Images are loaded into MATLAB by means of dedicated read functions.

For instance, a Heidelberg ``vol`` file can be read as:

.. code-block::

        [header, seg, bscan, fundus] = read_vol(file)

This function returns four arguments:

* ``header``: struct with metadata (image sdimensions, lateraliry, ...)
* ``seg``: struct with retinal layer segmentation
* ``bscan``: a 3D image with B-scans
* ``fundus``: a 2D image with a fundus image

Other file formats can be read similarly:

.. code-block:: matlab

        [header, seg, bscan, fundus] = read_e2e(file)

A single ``E2E`` file can contain more than one OCT volumes (usually two: one per eye).
For this reason, the returned objects are cell lists each related to a single OCT volume.

Topcon ``fda`` files contain metadata, OCT, segmentation and a color fundus image. To read them:

.. code-block:: matlab

        [header, seg, bscan, fundus] = read_fda(file)

Cirrus ``img`` files only contain 1 OCT volume and no segmentation or metadata.
The ``read_img`` function reads the B-scans and tries to infer the metadata from image size and file name.

.. code-block:: matlab

        [header, bscan] = read_img(file)

The fundus image acquired along with the img is stored in a separated ``bin`` file that can be read by:

.. code-block:: matlab

        fundus = read_binfile)

2. Explore metadata
^^^^^^^^^^^^^^^^^^^
After reading a file, it is a good idea to look into the metadata stored in the header to check that everything looks good.

We can do that with the following command:

.. code-block:: matlab

        disp(header)

Header structure is an attempt at harmonizing data from different scanners. 

Fields should be equivalent across scanners. However, due to the difficulty of fully reading propietary file formats, not all the fields will be present for every format.

We can also check which segmented layer boundaries are present in the file:

.. code-block:: matlab

        boundaries = fields(seg);
        disp(boundaries')

3. Visualize the images
^^^^^^^^^^^^^^^^^^^^^^^
Quality assurance is a key step in any image analysis pipeline. In OCT, poor contrast, artifacts, and segmentation errors are not uncommon.
Therefore, it is crucial to inspect the data to detect those problems.

Here we plot the middle B-scan with its segmentation:

.. code-block:: matlab

        f = figure;
        idx_bscan = 13;
        imshow(bscan(:,:,idx_bscan)); hold on;
        for j=1:length(boundaries)
            plot(seg.(boundaries{j})(idx_bscan, :))
        end


.. image:: images/bscan.png
    :align: center
    :width: 400
    :height: 370

Often a fundus image is acquired alongside the OCT volume. Depending on the scanner this image is different:

* Heidelberg: a graysacle image stored in ``vol``, ``e2e`` files.
* Cirrus: a grayscale image stored in a separated ``bin`` file
* Topcon: an RGB color image stored inside ``fda`` files.

To visualize the fundus image we can just:

.. code-block:: matlab
        
        f = figure;
        imagesc(fundus); colormap(gray); axis("off");

.. image:: images/fundus.png
    :align: center
    :width: 400
    :height: 400

4. Generate a report
^^^^^^^^^^^^^^^^^^^^
Doing all the previous steps separately is not always convenient.

To simplify data inspection, we can automatically create a summary report using the ``generate_report()`` function.
The function allows us to pass fundus,bscan and segmentation data and builds a tiled plot that can be stored in memory.

A full report includes:

* *Fundus image*  
* *Reflectance map*: a 2D map obtained by averaging A-scan intensities in depth. Useful to detect noisy B-scans and shades.
* *Thickness maps*: 2D maps generated by computing layer thicknesses from segmentation data. Useful to detect segmentation errors or lessions.
* *B-scans*: several slices with overlapping segmentation.

As an example, let's generate a full report for the loaded image.

As arguments we need to provide image data and the layers for which we want to plot thickness maps

.. code-block:: matlab

        layers = {'TRT','RNFL','GCIPL'};
        generate_report(bscan, seg, fundus, layers);

.. image:: images/report.png
    :align: center

Finally, to store it as a ``png`` withoug opening a window:

.. code-block:: matlab

        generate_report(bscan, seg, fundus, layers, 'file_name', 'report.png', 'visible', 'off');




**Feature extraction**
----------------------

.. _tutorial_feature_extraction:

This tutorial shows how to build a feature extraction pipeline step by step


0. Preliminaries
^^^^^^^^^^^^^^^^

The computation of retinal features requires the segmentation or retinal layer boundaries.

Here, there are two possiblities:

* The images were already segmented and the segmentation is stored in vol, e2e, or fda files
* We have raw images without segmentation

We can check that by reading the images and trying to visualize the segmentation as in the previous tutorial.

In the segmentation is not available, we recommend using `OCT Explorer <https://iibi.uiowa.edu/oct-reference>`_ software to do so. OCT Explorer allows you to load your images and save segmentation in a ``xml`` file that can be read by RETIMAT.

After we have the segmentation data we need to decide which features we want to compute. There are four categories of features:

* Thickness: the most common.
* Foveal pit morphology
* Reflectance: image quality
* Texture analysis

1. Load the files
^^^^^^^^^^^^^^^^^
This step parses the information inside the file and computes the ``X``, ``Y`` coordinates of each A-scan.

To analyse retinal thickness we only need the header and the segmentation data.

.. code-block:: matlab

        [header, seg, ~, ~] = read_vol(file, 'get_coordinates');

        X = header.X_oct;
        Y = header.Y_oct;

2. Preprocessing
^^^^^^^^^^^^^^^^

3. Feature computation
^^^^^^^^^^^^^^^^^^^^^^

4. Building the pipeline
^^^^^^^^^^^^^^^^^^^^^^^^

