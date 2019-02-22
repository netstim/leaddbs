#### The Advanced Examples
This examples require multiple datasets, e.g. go registered MRI or atlases. We currently cannot provide a complete dataset for reasons of public sharing data sharing permission. 
 * **EXAMPLE_1.m** - Basic PaCER call and electrode plot. Start here!
    * **EXAMPLE_1_1.m** - Continues EXAMPLE_1 by adding an **MPR view** of the CT image and demonstrating some **plot customisations**
    * **EXAMPLE_1_2.m** - Continues EXAMPLE_1 by adding **transformation of the electrode object to native T1 space** and loading MPR view of the native T1 image. This transformation is achieved by applying a FSL FLIRT transformation (from CT to T1 space) to the electrode model using the methods provided by the PaCER electrode objects.
 * **EXAMPLE_2.m** - Demonstrates the electrodes reconstruction from **longitudinal datasets** (immediate post-op, later post-op) as well as **visualisation of co-registered MRI images** (T1) as well as **visualisation of STN segmentations**. (cf. Use-Case A, Use-Case C).
 	* **EXAMPLE_2_1.m** - Continues EXAMPLE_2 by adding a **simplified VTA model** (Mädler/Coenen) (cf. Use-Case C)
 * **EXAMPLE_3.m** - Demonstrates plan/outcome comparisons by loading and displaying Medtronic(R) Framelink (TM) stereotactic plan with a PaCER electrode reconstructions (cf. Use-Case B).
 * **EXAMPLE_4.m** - Demonstrates PaCER operating in atlas space. Post OP CT and T1 linearly pre-registered to the template and electrodes plotted together with subcortical atlas segmentations. (cf. Use-Case D)
 * **EXAMPLE_5.m** - Demonstrates electrode reconstruction in **native** space with subsequent **transformation of the electrode objects to T1 space** (see example_1_2). The subcortical atlas from example 4 is now nonlinearly transformed and plotted with electrodes in T1 space. (cf. Use-Case D)
