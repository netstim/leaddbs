# from GUI import Ui_MainWindow
# from AppUI import MainWindow
from functionalities import Functionalities
from pop_up_control.cpe_active import PopUpCPEActive
from pop_up_control.external_neuron_arrray import PopUpExternalNeuronArray
from pop_up_control.full_field_ifft import PopUpFullFieldIFFT
from pop_up_control.mesh_refinement import PopUpMeshRefinement


class Dictionary(Functionalities):
    def __init__(self, win):
        self.main_win = win.main_win
        self.ui = win.ui

        # # call classes
        # self.externalNeuronArray = PopUpExternalNeuronArray(self.main_win)
        # self.meshRefinement = PopUpMeshRefinement(self.main_win)
        # self.cpeActive = PopUpCPEActive(self.main_win)
        # self.fullFieldIFFT = PopUpFullFieldIFFT(self.main_win)

    def output_dict(self):
        output_dict = {
            'Segm_MRI_processed': self.corrector(1, self.ui.checkBox_Voxel_orr_MRI.checkState()),
            'DTI_processed': self.corrector(0, self.ui.checkBox_Voxel_orr_DTI.checkState()),
            'Init_neuron_model_ready': self.corrector(1, self.ui.checkBox_Init_Neuron_Model_Ready.checkState()),
            'Init_mesh_ready': self.corrector(1, self.ui.checkBox_Init_Mesh_Ready.checkState()),
            'Adjusted_neuron_model_ready': self.corrector(1, self.ui.checkBox_Adjusted_Neuron_Model_Ready.checkState()),
            'CSF_mesh_ready': self.corrector(1, self.ui.checkBox_CSF_Mesh_Ready.checkState()),
            'Adapted_mesh_ready': self.corrector(0, self.ui.checkBox_Adapted_Mesh_Ready.checkState()),
            'Signal_generated': self.corrector(0, self.ui.checkBox_Signal_Generation_Ready.checkState()),
            'Parallel_comp_ready': self.corrector(0, self.ui.checkBox_Parallel_Computing_Ready.checkState()),
            'Parallel_comp_interrupted': self.corrector(0,
                                                        self.ui.checkBox_Parallel_Computing_Interrupted.checkState()),
            'IFFT_ready': self.corrector(0, self.ui.checkBox_IFFT_Ready.checkState()),
            'MRI_data_name': self.set_default_text('', self.ui.lineEdit_MRI_Data_Name.text()),
            'MRI_in_m': self.corrector(0, self.ui.checkBox_MRI_m.checkState()),
            'DTI_data_name': self.set_default_text('', self.ui.lineEdit_DTI_Data_Name.text()),
            'DTI_in_m': self.corrector(0, self.ui.checkBox_DTI_m.checkState()),
            'CSF_index': self.set_default_values(0, self.ui.doubleSpinBox_CSF_Index.value()),
            'WM_index': self.set_default_values(0, self.ui.doubleSpinBox_WM_Index.value()),
            'GM_index': self.set_default_values(0, self.ui.doubleSpinBox_GM_Index.value()),
            'default_material': self.set_default_values(2, self.ui.comboBox_Default_Material.currentIndex() + 1),
            'Electrode_type': self.set_default_text('SNEX100', self.ui.comboBox_Electrode_Type.currentText()),
            'Brain_shape_name': self.set_default_text('0', self.ui.lineEdit_Brain_Shape.text()),
            'Aprox_geometry_center': self.set_default_values(0, self.ui.lineEdit_Approx_Geometry_Center.text()),
            'Approximating_Dimensions': self.set_default_values(0, self.ui.lineEdit_Approximating_Dimensions.text()),
            'Implantation_coordinate_X': self.set_default_values(0,
                                                                 self.ui.doubleSpinBox_Implantation_Coordinate_X.value()),
            'Implantation_coordinate_Y': self.set_default_values(0,
                                                                 self.ui.doubleSpinBox_Implantation_Coordinate_Y.value()),
            'Implantation_coordinate_Z': self.set_default_values(0,
                                                                 self.ui.doubleSpinBox_Implantation_Coordinate_Z.value()),
            'Second_coordinate_X': self.set_default_values(0, self.ui.doubleSpinBox_2nd_Point_On_Lead_X.value()),
            'Second_coordinate_Y': self.set_default_values(0, self.ui.doubleSpinBox_2nd_Point_On_Lead_Y.value()),
            'Second_coordinate_Z': self.set_default_values(0, self.ui.doubleSpinBox_2nd_Point_On_Lead_Z.value()),
            'Rotation_Z': self.set_default_values(0, self.ui.doubleSpinBox_Turn_Around_Lead_Axis.value()),
            'encap_thickness': self.set_default_values(0, self.ui.doubleSpinBox_Encapsulation_Thickness.value()),
            'encap_tissue_type': self.ui.comboBox_Encapsulation_Tissue_Type.currentIndex() + 1,
            'encap_scaling_cond': self.ui.doubleSpinBox_Conductivity_Scaling.value(),
            'encap_scaling_perm': self.ui.doubleSpinBox_Permittivity_Scaling.value(),
            'pattern_model_name': self.set_default_text('0', self.ui.lineEdit_Pattern_Model_Name.text()),
            'Axon_Model_Type': self.set_default_text('McIntyre2002', self.ui.comboBox_Axon_Model_Type.currentText()),
            'diam_fib': self.set_default_values(0, self.ui.lineEdit_Fiber_Diameter.text()),
            'n_Ranvier': self.set_default_values(0, self.ui.lineEdit_N_Ranvier.text()),
            'v_init': self.set_default_values(0, self.ui.doubleSpinBox_V_Init.value()),
            'Neuron_model_array_prepared': self.corrector(0, self.ui.checkBox_Neuron_Model_Array_Prepared.checkState()),
            'Name_prepared_neuron_array': self.set_default_text('',
                                                                self.ui.lineEdit_Name_Prepared_Neuron_Array.text()),
            'EQS_core': self.set_default_text("1", self.ui.comboBox_Laplace_Formulation.currentText()),
            'Skip_mesh_refinement': self.corrector(0, self.ui.checkBox_Skip_Mesh_Refinement.checkState()),
            'el_order': self.set_default_values(0, self.ui.spinBox_El_Order.value()),
            'number_of_processors': self.set_default_values(0, self.ui.spinBox_Number_Of_Processors.value()),
            'FEniCS_MPI': self.corrector(0, self.ui.checkBox_FEniCS_MPI.checkState()),
            'current_control': self.corrector(0, self.ui.checkBox_Current_Control.checkState()),
            'Pulse_amp': self.set_default_values([0], self.ui.lineEdit_Phi_Vector.text()),
            'Solver_Type': self.set_default_text("MUMPS", self.ui.comboBox_Solver_Type.currentText()),
            'freq': self.set_default_values(130.0, self.ui.doubleSpinBox_Frequency.value()),
            'T': self.set_default_values(10.0e-5, self.ui.doubleSpinBox_T.value()),
            't_step': self.set_default_values(1 * 10e-7, self.ui.doubleSpinBox_T_Step.value()),
            'phi': self.set_default_values(0.0, self.ui.doubleSpinBox_Signal_Shift.value()),
            'Signal_type': self.set_default_text("Rectangle", self.ui.comboBox_Signal_Type.currentText()),
            'Ampl_scale': self.set_default_values(1.0, self.ui.doubleSpinBox_Amplitude_Scale.value()),
            'CPE_activ': self.corrector(0, self.ui.checkBox_CPE_Active.checkState()),
            'VTA_approx': self.corrector(0, self.ui.checkBox_Full_Field_IFFT.checkState()),
            'spectrum_trunc_method': self.set_default_text("No Truncation",
                                                           self.ui.comboBox_Spectrum_Truncation_Method.currentText()),

            'trunc_param': self.set_default_values(1, self.spectrum_truncation_method_output(
                'Octave Band Method', self.ui.comboBox_Spectrum_Truncation_Method,
                self.ui.lineEdit_Truncation_Parameter, self.ui.spinBox_Truncation_Parameter
            )),

            'Truncate_the_obtained_full_solution': self.corrector(0,
                                                                  self.ui.checkBox_Truncate_The_Obtained_Full_Solution.checkState()),
            #'Show_paraview_screenshots': self.corrector(0, self.ui.checkBox_Show_Paraview_Screenshots.checkState()),
            'external_grounding': self.corrector(0, self.ui.checkBox_external_grounding.checkState())
        }
        return output_dict
