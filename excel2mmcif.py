# This script reads in a excel file of a fixed format and uses the python-ihm
# implementation to generate the objects and write the mmcif file in the end
# Christian Hanke 11.03.2019
# christian.hanke@hhu.de
# version 1.06

# Note: In case of non-mandatory parameters, it should be checked whether the
# column is present at all or whether the respective cell in the excel sheet
# is empty (using pandas.isnull).

# Note: Currently, the script requires to be used with Python 3 in order to
# ensure proper output in the mmcif files of the unicode strings.
# Note: The script requires the python ihm, pandas, and openpyxl

import pandas
import argparse
import ihm
import ihm.flr
import ihm.dumper
import ihm.dataset
import ihm.location
import ihm.startmodel
import ihm.analysis
import ihm.protocol
import ihm.model
import ihm.source
import ihm.reader

BEVERBOSE = True
DEBUG = False


class DummyClass():
    def __init__(self, value):
        self._id = value


def occurs_in_list(curobject, list):
    for entry in list:
        if curobject.__dict__ == entry.__dict__:
            return True
    return False


def get_index_from_list(curobject, list):
    for i in range(len(list)):
        if curobject.__dict__ == list[i].__dict__:
            return i
    return -1

def get_resatom_from_list(curobject, list):
    """ Get the residue or atom from a list by checking for identical entries
    in the list. This is not achieved by simply comparing, since Atom and
    Residue use __slots__.
    In case of an Atom, also the respective Residue has to be checked for
    equality.
    :param curobject: Atom or Residue
    :param list: The resatom list
    :return: The object that has all the same attribute values as curobject
        or None.
    """
    for i in range(len(list)):
        # check whether all attribute values of the objects are identical
        isidentical = True
        # Check whether the objects have the same attrributes
        if curobject.__slots__ == list[i].__slots__:
            for attr in curobject.__slots__:
                curobjectattrvalue = getattr(curobject, attr)
                listobjectattrvalue = getattr(list[i], attr)
                # in case there is an additional object with slots in the
                # curobject we also need to check the slots of this one
                if hasattr(curobjectattrvalue, '__slots__'):
                    if curobjectattrvalue.__slots__ == \
                            listobjectattrvalue.__slots__:
                        for innerattr in curobjectattrvalue.__slots__:
                            # Check whether innerattr is included in an object
                            if hasattr(curobjectattrvalue, innerattr):
                                # and in the other object
                                if hasattr(listobjectattrvalue, innerattr):
                                    # Check whether the attributes are the
                                    # same
                                    if not getattr(curobjectattrvalue,
                                                   innerattr) == \
                                           getattr(listobjectattrvalue,
                                                   innerattr):
                                        isidentical = False
                                # If only one of the objects has the
                                # attributes, they are not identical
                                else:
                                    isidentical = False
                # otherwise, we just compare them
                else:
                    if not (curobjectattrvalue == listobjectattrvalue):
                        isidentical = False
            if isidentical:
                return list[i]
    return None


class Excel2mmcifConverter:
    def __init__(self,
                 excel_filename,
                 cifout_filename,
                 atom_site_filename,
                 be_verbose=True):
        self.verbose=be_verbose
        self.excel_filename = excel_filename
        self.cifout_filename = cifout_filename
        self.atom_site_filename = atom_site_filename

        # Create an IHM system
        self.system = ihm.System()
        # Create an excel file object
        self.xls_file = pandas.ExcelFile(self.excel_filename)
        # Prepare the lists to store data
        self.initialize_lists()

    def run(self):
        """Read and process the data from the excel sheet, and write the cif file"""
        # Fill all the data
        self.fill()
        # Write the file
        self.write_cif_file()

    def initialize_lists(self):
        """Prepare lists and dictionaries to store data"""
        # Collection of information to be stored
        self.list_resatoms = []

        self.list_citations = []
        self.list_citation_ids = []
        self.list_ihm_softwares = []
        self.list_ihm_software_ids = []
        self.list_locations = []
        # Datasets
        self.list_datasets = []
        self.list_dataset_ids = []
        self.list_dataset_groups = []
        self.list_dataset_group_ids = []
        # Dataset subgroups correspond to the entries in the Dataset list
        # (i.e. there might be a subgroup for Single-molecule FRET data and
        # another one for De Novo models)
        self.list_dataset_subgroup_by_dataset_list_id = []
        self.list_dataset_subgroup_by_dataset_list_id_ids = []
        # External reference ids for connection to the external files
        self.list_dataset_external_reference_ids = []
        # connect datasets with external files; store datasets and the
        # respective external_file ids
        self.list_connection_dataset_external_file_datasets = []
        self.list_connection_dataset_external_file_ids = []
        # repositories
        self.list_repositories = []
        self.list_respository_ids = []
        # External files
        self.list_external_files = []
        self.list_external_files_locations = []
        self.list_external_files_ids = []
        self.tmp_dict_for_dataset_groups = {}
        # Entity (IHM) and Entity assembly (FLR)
        self.list_ihm_entities = []
        self.list_ihm_entity_ids = []
        self.list_entity_assemblies = []
        self.list_entity_assembly_ids = []
        # Comparative starting model
        self.list_starting_model_templates = []
        self.tmp_dict_templates_for_starting_models = {}
        # Instance (AsymUnit)
        self.list_asym_units = []
        self.list_asym_units_ids = []
        self.list_asym_unit_ranges = []
        self.list_ihm_model_representations = []
        self.list_ihm_model_representation_ids = []
        self.list_ihm_starting_model_details = []
        self.list_structure_assemblies = []
        self.list_structure_assembly_ids = []
        self.tmp_dict_for_structure_assembly = {}
        self.tmp_dict_for_structure_assembly_names = {}
        self.tmp_dict_for_structure_assembly_description = {}
        self.tmp_dict_for_model_representations = {}
        # Modeling protocol and modeling post process
        self.list_ihm_modeling_protocol_analysis_steps = []
        self.list_ihm_modeling_protocol_analysis_step_ids = []
        self.list_ihm_modeling_protocols = []
        self.list_ihm_modeling_protocols_ids = []
        self.list_ihm_modeling_protocol_analyses = []
        self.list_ihm_modeling_protocol_analyses_ids = []
        self.list_ihm_modeling_protocol_modeling_steps = []
        self.list_ihm_modeling_protocol_modeling_step_ids = []
        self.dict_ihm_modeling_protocol_software_ids = {}
        # Multi-state modeling
        self.list_empty_states = []
        self.list_empty_state_ids = []
        # Models
        self.list_models = []
        self.list_models_ids = []
        # Models - Model groups
        self.list_model_groups = []
        self.list_model_group_ids = []
        # Models - models that belong to one state
        self.list_models_states = []
        self.list_models_state_ids = []
        self.list_models_state_groups = []
        self.list_models_state_group_ids = []
        # Ensembles
        self.list_ensembles = []
        # FLR
        # FLR - FPS global parameters
        self.list_flr_fps_global_parameters = []
        # A dictionary to identify the global parameters by the protocol id.
        # The dictionary, which uses the protocol id as keys contains additional
        # dictionaries that use the protocol steps as keys
        self.dict_flr_fps_global_parameters_by_protocol_id = {}
        # FLR - FPS MPP
        self.dict_flr_fps_mpp_groups = {}
        # The following lists are for both, reference measurements and FLR
        # FLR - sample conditions
        self.list_sample_conditions = []
        # FLR - instruments
        self.list_instruments = []
        # FLR - instrument settings
        self.list_inst_settings = []
        # FLR - experimental conditions
        self.list_exp_conditions = []
        # FLR - samples
        self.list_samples = []
        self.list_sample_ids = []
        # FLR - experiments
        self.list_experiments = []
        # FLR - chemical descriptors
        self.list_chemical_descriptors = []
        # FLR - chemical components
        self.list_chemical_components = []
        # FLR - probes
        self.list_probe_donors = []
        self.list_probe_acceptors = []
        # FLR - modified and mutated residues
        self.list_modified_residues = []
        self.list_mutated_residues = []
        # Poly_probe_positions
        self.list_poly_probe_positions = []
        self.list_residue_objects = []
        # Sample_probe_details
        self.list_sample_probe_details = []
        # Probe_descriptor_conjugate
        self.list_probe_conjugate_descriptors = []
        # Poly_probe_conjugate
        self.list_poly_probe_conjugates = []
        # Reference measurement lifetimes and fractions
        self.list_ref_measurement_lifetimes = []
        # Reference measurements
        self.list_ref_measurements = []
        # Reference measurement groups
        self.list_ref_measurement_groups = []
        self.list_ref_measurement_group_ids = []
        self.list_of_object_indices_refmeas = []
        # FLR
        self.list_of_object_indices = []
        self.list_fret_forster_radius = []
        self.list_fret_calibration_parameters = []
        self.list_lifetime_fit_models = []
        self.list_fret_analysis = []
        self.list_peak_assignments = []
        self.list_fret_distance_restraint_groups = []
        self.list_fret_distance_restraints = []
        self.list_fret_distance_restraint_ids = []
        self.list_FPS_AV_parameters = []
        self.list_FPS_modeling_by_AV = []
        self.list_FPS_AV_modeling = []
        self.list_FPS_mean_probe_positions = []
        self.list_FPS_modeling_by_MPP = []
        self.list_FPS_MPP_modeling = []
        self.list_flr_model_quality = []
        self.list_flr_model_distances = []

    def fill(self):
        self.add_general()
        self.add_citation()
        self.add_ihm_software()
        self.add_ihm_dataset_and_external_files()
        self.add_ihm_entity()
        self.add_ihm_comparative_starting_model()
        self.add_ihm_instance_and_AsymUnits()
        self.add_ihm_modeling_protocol()
        self.add_ihm_multi_state_modeling_and_models()
        self.add_ihm_ensemble()
        self.add_flr_fps_global_parameters()
        self.add_flr_fps_mpp()
        self.add_flr_reference_measurements()
        self.add_flr()
        self.add_flr_model_quality()
        self.add_flr_model_distances()
        self.create_flrdata()
        self.handle_atom_site_file()

    def add_general(self, sheet_name='General', skiprows=3, header=0):
        """Read the general information from the excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'General\' ...")
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_Struct_id'])
        for i in range(nr_entries):
            cur_general_title = \
                None if ('IHM_Struct_title' not in xls_data.keys()
                         or pandas.isnull(xls_data['IHM_Struct_title'][i])) \
                    else xls_data['IHM_Struct_title'][i]
            self.system.title = cur_general_title

    def add_citation(self, sheet_name='Citation', skiprows=3, header=0):
        """Read the information on the Citation"""
        if self.verbose:
            print(" ... Processing tab \'Citation\' ...")
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_Citation_id'])

        for i in range(nr_entries):
            cur_citation_id = xls_data['IHM_Citation_id'][i]
            cur_citation_title = xls_data['IHM_Citation_Title'][i]
            cur_citation_journal_abbreviation = \
                None if ('IHM_Citation_Journal_abbreviation' not in
                         xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Citation_Journal_abbreviation'][i])) \
                    else xls_data['IHM_Citation_Journal_abbreviation'][i]
            cur_citation_journal_volume = \
                None if ('IHM_Citation_Journal_volume' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Citation_Journal_volume'][i])) \
                    else xls_data['IHM_Citation_Journal_volume'][i]
            cur_citation_first_page = \
                None if ('IHM_Citation_First_page' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Citation_First_page'][i])) \
                    else xls_data['IHM_Citation_First_page'][i]
            cur_citation_last_page = \
                None if ('IHM_Citation_Last_page' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Citation_Last_page'][i])) \
                    else xls_data['IHM_Citation_Last_page'][i]
            cur_citation_year = \
                None if ('IHM_Citation_Year' not in xls_data.keys()
                         or pandas.isnull(xls_data['IHM_Citation_Year'][i])) \
                    else xls_data['IHM_Citation_Year'][i]
            cur_citation_pubmed_id = \
                None if ('IHM_Citation_Pubmed_id' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Citation_Pubmed_id'][i])) \
                    else xls_data['IHM_Citation_Pubmed_id'][i]
            cur_citation_doi = \
                None if ('IHM_Citation_DOI' not in xls_data.keys()
                         or pandas.isnull(xls_data['IHM_Citation_DOI'][i])) \
                    else xls_data['IHM_Citation_DOI'][i]
            cur_citation_authors = xls_data['IHM_Citation_Authors'][i]

            cur_citation_authors_list = cur_citation_authors.split(';')
            # Create the citation object
            cur_citation = ihm.Citation(
                pmid=cur_citation_pubmed_id,
                title=cur_citation_title,
                journal=cur_citation_journal_abbreviation,
                volume=cur_citation_journal_volume,
                page_range=(cur_citation_first_page, cur_citation_last_page),
                year=cur_citation_year,
                authors=cur_citation_authors_list,
                doi=cur_citation_doi)

            if not occurs_in_list(cur_citation, self.list_citations):
                self.list_citations.append(cur_citation)
                self.list_citation_ids.append(cur_citation_id)

        # add all citation objects to the system
        for entry in self.list_citations:
            self.system.citations.append(entry)

    def add_ihm_software(self, sheet_name='Software', skiprows=3, header=0):
        """Read information on the Software from the excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'Software\' ...")
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_Software_id'])

        for i in range(nr_entries):
            cur_ihm_software_id = xls_data['IHM_Software_id'][i]
            cur_ihm_software_name = xls_data['IHM_Software_name'][i]
            cur_ihm_software_classification = \
                xls_data['IHM_Software_classification'][i]
            cur_ihm_software_description = \
                None if ('IHM_Software_description' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Software_description'][i])) \
                    else xls_data['IHM_Software_description'][i]
            cur_ihm_software_location = \
                None if ('IHM_Software_location' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Software_location'][i]))\
                    else xls_data['IHM_Software_location'][i]
            cur_ihm_software_type = xls_data['IHM_Software_type'][i]
            cur_ihm_software_version = \
                None if ('IHM_Software_version' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Software_version'][i])) \
                    else xls_data['IHM_Software_version'][i]
            cur_ihm_software = ihm.Software(
                name=cur_ihm_software_name,
                classification=cur_ihm_software_classification,
                description=cur_ihm_software_description,
                location=cur_ihm_software_location,
                type=cur_ihm_software_type,
                version=cur_ihm_software_version)
            if cur_ihm_software not in self.list_ihm_softwares:
                self.list_ihm_softwares.append(cur_ihm_software)
                self.list_ihm_software_ids.append(cur_ihm_software_id)

        for entry in self.list_ihm_softwares:
            self.system.software.append(entry)

    def add_ihm_dataset_and_external_files(self,
                                           sheet_name_dataset='Dataset',
                                           skiprows_dataset=3,
                                           header_dataset=0,
                                           sheet_name_external_files=
                                           'External_files',
                                           skiprows_external_files=3,
                                           header_external_files=0):
        """Read on the Datasets and external file information from the
        excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'Dataset\' ...")

        # Note: Many of these temporary lists could possibly be avoided when
        # using more classes. => Future TODO
        # storage for datasets that belong to one group and create a
        # dataset_group later
        tmp_dict_for_dataset_groups = {}
        tmp_info_for_dataset_groups = {}
        tmp_dict_for_dataset_subgroups = {}
        # store the type of the dataset for each external reference id. This
        # will be used for the external files to create the datasets
        tmp_dict_for_external_reference_store_dataset_type = {}
        tmp_dict_for_external_reference_store_dataset_group = {}
        tmp_dict_for_external_reference_store_dataset_details = {}
        tmp_dict_for_external_reference_store_repository = {}
        tmp_dict_for_external_reference_store_dataset_list_id = {}
        # Read the data
        xls_data_d = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name_dataset,
                                     skiprows=skiprows_dataset,
                                     header=header_dataset)
        nr_entries_d = len(xls_data_d['IHM_Dataset_Dataset_list_id'])
        for i in range(nr_entries_d):
            cur_dataset_list_id = xls_data_d['IHM_Dataset_Dataset_list_id'][i]
            # required to create the dataset_group
            cur_dataset_group = xls_data_d['IHM_Dataset_Dataset_group'][i]
            cur_dataset_group_name = \
                None if ('IHM_Dataset_Dataset_group_name' not in
                         xls_data_d.keys() or pandas.isnull(
                            xls_data_d['IHM_Dataset_Dataset_group_name'][i])) \
                    else xls_data_d['IHM_Dataset_Dataset_group_name'][i]
            cur_dataset_group_details = \
                None if (
                        'IHM_Dataset_Dataset_group_details' not in
                        xls_data_d.keys() or pandas.isnull(
                    xls_data_d['IHM_Dataset_Dataset_group_details'][i])) \
                    else xls_data_d['IHM_Dataset_Dataset_group_details'][i]
            # required to create the dataset of the correct type
            cur_dataset_data_type = xls_data_d['IHM_Dataset_Data_type'][i]
            # required for the creation of the database location
            cur_dataset_DB_flag = xls_data_d['IHM_Dataset_DB_flag'][i] in \
                                  ['Yes', 'YES', 'yes']
            cur_dataset_DB_name = xls_data_d['IHM_Dataset_DB_name'][i]
            cur_dataset_DB_accession_code = \
                xls_data_d['IHM_Dataset_DB_accession_code'][i]
            cur_dataset_DB_version = \
                None if ('IHM_Dataset_DB_version' not in xls_data_d.keys()
                         or pandas.isnull(
                            xls_data_d['IHM_Dataset_DB_version'][i])) \
                    else xls_data_d['IHM_Dataset_DB_version'][i]
            # required for the connection to the external files
            cur_dataset_external_reference_id = \
                xls_data_d['IHM_Dataset_External_reference_id'][i]
            # required for the creation of the repository
            cur_dataset_reference_provider = \
                xls_data_d['IHM_Dataset_Reference_provider'][i]
            cur_dataset_reference_type = \
                xls_data_d['IHM_Dataset_Reference_type'][i]
            cur_dataset_reference = xls_data_d['IHM_Dataset_Reference'][i]
            cur_dataset_refers_to = xls_data_d['IHM_Dataset_Refers_to'][i]
            cur_dataset_associated_url = \
                xls_data_d['IHM_Dataset_Associated_url'][i]
            cur_dataset_top_directory = \
                None if ('IHM_Dataset_top_directory' not in xls_data_d.keys()
                         or pandas.isnull(
                            xls_data_d['IHM_Dataset_top_directory'][i])) \
                    else xls_data_d['IHM_Dataset_top_directory'][i]
            cur_dataset_details = \
                None if ('IHM_Dataset_Details' not in xls_data_d.keys()
                         or pandas.isnull(
                            xls_data_d['IHM_Dataset_Details'][i])) \
                    else xls_data_d['IHM_Dataset_Details'][i]

            # if the dataset is a database entry
            if cur_dataset_DB_flag:
                cur_dataset = None
                cur_location = None
                # PDB
                if cur_dataset_DB_name == 'PDB':
                    cur_location = ihm.location.PDBLocation(
                        db_code=cur_dataset_DB_accession_code,
                        version=cur_dataset_DB_version,
                        details=cur_dataset_details)
                    if cur_location is not None and \
                            cur_location not in self.list_locations:
                        self.list_locations.append(cur_location)
                    cur_dataset = ihm.dataset.PDBDataset(cur_location)

                # SASBDB
                elif cur_dataset_DB_name == 'SASBDB':
                    cur_location = ihm.location.SASBDBLocation(
                        db_code=cur_dataset_DB_accession_code,
                        version=cur_dataset_DB_version,
                        details=cur_dataset_details)
                    if cur_location not in self.list_locations:
                        self.list_locations.append(cur_location)
                    cur_dataset = ihm.dataset.SASDataset(cur_location)

                # BMRB
                elif cur_dataset_DB_name == 'BMRB':
                    cur_location = ihm.location.BMRBLocation(
                        db_code=cur_dataset_DB_accession_code,
                        version=cur_dataset_DB_version,
                        details=cur_dataset_details)
                    if cur_location not in self.list_locations:
                        self.list_locations.append(cur_location)
                    cur_dataset = ihm.dataset.NMRDataset(cur_location)

                # TODO: handle other databases
                else:
                    print('NOTE! Database name %s not handled.' % (
                        cur_dataset_DB_name))

                if cur_dataset is not None:
                    if cur_dataset not in self.list_datasets:
                        self.list_datasets.append(cur_dataset)
                        self.list_dataset_ids.append(cur_dataset_list_id)
                        self.list_dataset_external_reference_ids.append(
                            cur_dataset_external_reference_id)
                    if cur_ihm_dataset_group not in \
                            tmp_dict_for_dataset_groups.keys():
                        tmp_dict_for_dataset_groups[cur_dataset_group] = []
                        tmp_dict_for_dataset_groups[cur_dataset_group].append(
                            cur_dataset)
                        tmp_dict_for_dataset_groups[cur_dataset_group] = {}
                        tmp_info_for_dataset_groups[cur_dataset_group]['name'] \
                            = cur_ihm_dataset_group_name
                        tmp_info_for_dataset_groups[
                            cur_dataset_group]['details'] = \
                            cur_dataset_group_details
                    else:
                        tmp_dict_for_dataset_groups[cur_dataset_group].append(
                            cur_dataset)

            # otherwise it is stored in a repository
            else:
                cur_dataset = None
                if cur_dataset_reference_type == 'DOI':
                    ## !!! TODO: root
                    cur_repo = ihm.location.Repository(
                        doi=cur_dataset_reference,
                        root='.',
                        url=cur_dataset_associated_url,
                        top_directory=cur_dataset_top_directory)
                    if not occurs_in_list(cur_repo, self.list_repositories):
                        self.list_repositories.append(cur_repo)
                        tmp_dict_for_external_reference_store_dataset_list_id[
                            cur_dataset_external_reference_id] = cur_dataset_list_id
                        tmp_dict_for_external_reference_store_dataset_type[
                            cur_dataset_external_reference_id] = cur_dataset_data_type
                        tmp_dict_for_external_reference_store_dataset_group[
                            cur_dataset_external_reference_id] = cur_dataset_group
                        tmp_dict_for_external_reference_store_repository[
                            cur_dataset_external_reference_id] = cur_repo
                        tmp_dict_for_external_reference_store_dataset_details[
                            cur_dataset_external_reference_id] = cur_dataset_details

                    # Add dataset group id to the list of dataset groups;
                    # Still needs to be filled.
                    if cur_dataset_group not in \
                            tmp_dict_for_dataset_groups.keys():
                        tmp_dict_for_dataset_groups[cur_dataset_group] = []
                        tmp_info_for_dataset_groups[cur_dataset_group] = {}
                        tmp_info_for_dataset_groups[cur_dataset_group]['name'] = \
                            cur_dataset_group_name
                        tmp_info_for_dataset_groups[cur_dataset_group]['details'] = \
                            cur_dataset_group_details

        # Now read the external file information
        if self.verbose:
            print(" ... Processing tab \'External_files\' ... ")
        xls_data_e = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name_external_files,
                                     skiprows=skiprows_external_files,
                                     header=header_external_files)
        nr_entries_e = len(xls_data_e['IHM_External_file_Ordinal'])
        for i in range(nr_entries_e):
            cur_external_files_ordinal = \
                xls_data_e['IHM_External_file_Ordinal'][i]
            cur_external_files_reference_id = \
                xls_data_e['IHM_External_file_Reference_id'][i]
            cur_external_files_file_path = \
                xls_data_e['IHM_External_file_File_path'][i]
            cur_external_files_file_format = \
                xls_data_e['IHM_External_file_File_format'][i]
            cur_external_files_content_type = \
                xls_data_e['IHM_External_file_Content_type'][i]
            cur_external_files_file_size = \
                xls_data_e['IHM_External_file_File_size'][i]
            cur_external_files_details = \
                None if ('IHM_External_file_Details' not in xls_data_e.keys()
                         or pandas.isnull(
                            xls_data_e['IHM_External_file_Details'][i])) \
                    else xls_data_e['IHM_External_file_Details'][i]

            # from the dataset tab
            cur_external_files_dataset_type = \
                tmp_dict_for_external_reference_store_dataset_type[
                    cur_external_files_reference_id]
            cur_external_files_dataset_group = \
                tmp_dict_for_external_reference_store_dataset_group[
                    cur_external_files_reference_id]
            cur_external_files_dataset_details = \
                tmp_dict_for_external_reference_store_dataset_details[
                    cur_external_files_reference_id]
            cur_external_files_repository = \
                tmp_dict_for_external_reference_store_repository[
                    cur_external_files_reference_id]
            cur_external_files_dataset_list_id = \
                tmp_dict_for_external_reference_store_dataset_list_id[
                    cur_external_files_reference_id]

            cur_location = None
            cur_dataset = None
            if cur_external_files_content_type == 'Input data or restraints':
                cur_location = \
                    ihm.location.InputFileLocation(
                        path=cur_external_files_file_path,
                        repo=cur_external_files_repository,
                        details=cur_external_files_details)
            if cur_external_files_content_type == \
                    'Modeling or post-processing output':
                cur_location = \
                    ihm.location.OutputFileLocation(
                        path=cur_external_files_file_path,
                        repo=cur_external_files_repository,
                        details=cur_external_files_details)
            if cur_external_files_content_type == \
                    'Modeling workflow or script':
                cur_location = \
                    ihm.location.WorkflowFileLocation(
                        path=cur_external_files_file_path,
                        repo=cur_external_files_repository,
                        details=cur_external_files_details)
            if cur_external_files_content_type == 'Visualization script':
                cur_location = \
                    ihm.location.VisualizationFileLocation(
                        path=cur_external_files_file_path,
                        repo=cur_external_files_repository,
                        details=cur_external_files_details)
            if cur_external_files_content_type == 'Other':
                cur_location = \
                    ihm.location.FileLocation(
                        path=cur_external_files_file_path,
                        repo=cur_external_files_repository,
                        details=cur_external_files_details)
            if cur_location not in self.list_external_files_locations:
                self.list_external_files_locations.append(cur_location)
                self.list_external_files_ids.append(cur_external_files_ordinal)
                self.system.locations.append(cur_location)

            # create the dataset
            cur_dataset = None
            if cur_external_files_dataset_type == 'Single molecule FRET data':
                cur_dataset = ihm.dataset.FRETDataset(
                    location=cur_location,
                    details=cur_external_files_dataset_details)
            elif cur_external_files_dataset_type == 'Ensemble FRET data':
                cur_dataset = ihm.dataset.EnsembleFRETDataset(
                    location=cur_location,
                    details=cur_external_files_dataset_details)
            elif cur_external_files_dataset_type == 'De Novo model':
                cur_dataset = ihm.dataset.DeNovoModelDataset(
                    location=cur_location,
                    details=cur_external_files_dataset_details)
            elif cur_external_files_dataset_type == 'Integrative model':
                cur_dataset = ihm.dataset.IntegrativeModelDataset(
                    location=cur_location,
                    details=cur_external_files_dataset_details)
            elif cur_external_files_dataset_type == 'Other':
                cur_dataset = ihm.dataset.Dataset(
                    location=cur_location,
                    details=cur_external_files_dataset_details)
            else:
                cur_dataset = ihm.dataset.Dataset(
                    location=cur_location,
                    details=cur_external_files_dataset_details)
            # store the dataset in the list for dataset groups
            if cur_dataset is not None and not occurs_in_list(
                    cur_dataset,
                    self.list_datasets):
                self.list_datasets.append(cur_dataset)
                self.list_dataset_ids.append(cur_external_files_dataset_list_id)
                # Store the external_files_reference_id as well in order to be
                # able to use it for the Dataset groups
                self.list_dataset_external_reference_ids.append(
                    cur_external_files_reference_id)
                if cur_external_files_dataset_group not in \
                        tmp_dict_for_dataset_groups.keys():
                    tmp_dict_for_dataset_groups[
                        cur_external_files_dataset_group] = []
                tmp_dict_for_dataset_groups[
                    cur_external_files_dataset_group].append(cur_dataset)
            # otherwise, use the previously generated dataset
            else:
                cur_dataset = self.list_datasets[get_index_from_list(
                    cur_dataset,
                    self.list_datasets)]
            if cur_external_files_dataset_group not in \
                    tmp_dict_for_dataset_groups.keys():
                tmp_dict_for_dataset_groups[
                    cur_external_files_dataset_group] = []
            if cur_dataset is not None:
                tmp_dict_for_dataset_groups[
                    cur_external_files_dataset_group].append(cur_dataset)
            # store the current dataset for the subgroups that correspond to the
            # dataset_list_ids
            if cur_external_files_reference_id not in \
                    tmp_dict_for_dataset_subgroups.keys():
                tmp_dict_for_dataset_subgroups[
                    cur_external_files_reference_id] = []
            if cur_dataset is not None:
                if cur_dataset not in \
                        tmp_dict_for_dataset_subgroups[
                            cur_external_files_reference_id]:
                    tmp_dict_for_dataset_subgroups[
                        cur_external_files_reference_id].append(cur_dataset)

            # store the connection between datasets and external file ordinal ids
            if cur_dataset is not None:
                self.list_connection_dataset_external_file_datasets.append(
                    cur_dataset)
                self.list_connection_dataset_external_file_ids.append(
                    cur_external_files_ordinal)

        # Go through the Dataset entries again to collect all the external
        # references belonging to one of the dataset groups
        # This is important because an external file could be used in multiple
        # dataset groups
        # Only if external files are present at all.
        for i in range(nr_entries_d):
            cur_dataset_group = \
                xls_data_d['IHM_Dataset_Dataset_group'][i]
            cur_dataset_DB_flag = \
                xls_data_d['IHM_Dataset_DB_flag'][i] in \
                ['Yes', 'YES', 'yes']
            cur_dataset_external_reference_id = \
                xls_data_d['IHM_Dataset_External_reference_id'][i]
            # Only if the dataset entry has an external file and is not in a
            # database
            if not cur_dataset_DB_flag:
                # if the current dataset group is not in the list of dataset
                # groups yet
                if not cur_dataset_group in tmp_dict_for_dataset_groups.keys():
                    # add it
                    tmp_dict_for_dataset_groups[cur_dataset_group] = []
                # then add the data
                # check for each of the previously collected entries in the
                # subgroups (which were collected by external reference id)
                # whether they are already added
                for cur_entry in tmp_dict_for_dataset_subgroups[
                    cur_dataset_external_reference_id]:
                    # check whether the entry is already in the list for the
                    # dataset groups
                    if cur_entry not in \
                            tmp_dict_for_dataset_groups[cur_dataset_group]:
                        tmp_dict_for_dataset_groups[
                            cur_dataset_group].append(cur_entry)

        # create the dataset_group
        for groupkey in tmp_dict_for_dataset_groups.keys():
            cur_dataset_group = \
                ihm.dataset.DatasetGroup(
                    elements=tmp_dict_for_dataset_groups[groupkey],
                    name=tmp_info_for_dataset_groups[groupkey]['name'],
                    details=tmp_info_for_dataset_groups[groupkey]['details'])
            self.list_dataset_groups.append(cur_dataset_group)
            self.list_dataset_group_ids.append(groupkey)

        # create the subgroups (corresponding to the dataset_list_id)
        for groupkey in tmp_dict_for_dataset_subgroups.keys():
            cur_dataset_subgroup = \
                ihm.dataset.DatasetGroup(tmp_dict_for_dataset_subgroups[groupkey])
            self.list_dataset_subgroup_by_dataset_list_id.append(cur_dataset_subgroup)
            self.list_dataset_subgroup_by_dataset_list_id_ids.append(groupkey)
        # update the locations in the repositories
        self.system.update_locations_in_repositories(self.list_repositories)

    def add_ihm_entity(self, sheet_name='Entity', skiprows=3, header=0):
        """Read entity information from the excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'Entity and Entity Assembly (FLR)\' "
                  "...")
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_Entity_Ordinal'])

        # Entity (IHM) and Entity assembly (FLR)
        for i in range(nr_entries):
            cur_entity_molecular_entity = \
                xls_data['IHM_Entity_Molecular_entity'][i]
            cur_entity_type = xls_data['IHM_Entity_Entity_type'][i]
            cur_entity_source_method = xls_data['IHM_Entity_Source_method'][i]
            cur_entity_description = \
                None if ('IHM_Entity_Description' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Entity_Description'][i])) \
                    else xls_data['IHM_Entity_Description'][i]
            cur_entity_polymer_type = xls_data['IHM_Entity_Polymer_type'][i]
            cur_entity_polymer_one_letter_code = \
                xls_data['IHM_Entity_Polymer_one_letter_code'][i]
            cur_entity_polymer_one_letter_code_canonical = \
                xls_data['IHM_Entity_Polymer_one_letter_code_canonical'][i]
            cur_entity_nonpolymer_chem_comp_id = \
                xls_data['IHM_Entity_Nonpolymer_chem_comp_ID'][i]
            cur_entity_nonpolymer_chem_comp_name = \
                xls_data['IHM_Entity_Nonpolymer_chem_comp_name'][i]
            cur_entity_nonpolymer_chem_comp_formula = \
                xls_data['IHM_Entity_Nonpolymer_chem_comp_formula'][i]

            # Source
            ## !!! TODO: For natural and synthetic: extend for the details
            cur_source = None
            if cur_entity_source_method == 'genetically manipulated source':
                cur_source = ihm.source.Manipulated()
            elif cur_entity_source_method == 'natural source':
                cur_source = ihm.source.Natural(Details=None)
            elif cur_entity_source_method == 'synthetic source':
                cur_source = ihm.source.Synthetic(Details=None)
            else:
                cur_source = ihm.source.Source()

            cur_entity = None
            # If the entity is a polymer
            if cur_entity_type == 'polymer':
                cur_alphabet = {
                    'polypeptide(D)': ihm.DPeptideAlphabet(),
                    'polypeptide(L)': ihm.LPeptideAlphabet(),
                    'polyribonucleotide': ihm.RNAAlphabet(),
                    'polydeoxribonucleotide':
                        ihm.DNAAlphabet()}[cur_entity_polymer_type]
                try:
                    cur_entity = ihm.Entity(
                        sequence=cur_entity_polymer_one_letter_code,
                        alphabet=cur_alphabet,
                        description=cur_entity_description,
                        source=cur_source)
                except KeyError:
                    # Check whether there are modified residues in the sequence
                    # Note: Modified residues have to be denoted by (Xyz),
                    # i.e. Three letters in brackets.
                    if '(' in cur_entity_polymer_one_letter_code:
                        list_of_sequence = []
                        this_index_i = 0  ## non-canonical index
                        this_index_j = 0  ## canonical index
                        ## Go through the non-canonical sequence
                        while this_index_i < \
                                len(cur_entity_polymer_one_letter_code):
                            # If non-canonical and canonical sequence match,
                            # this residue is kept
                            if cur_entity_polymer_one_letter_code[this_index_i] == \
                                    cur_entity_polymer_one_letter_code_canonical[this_index_j]:
                                list_of_sequence.append(
                                    cur_entity_polymer_one_letter_code[this_index_i])
                                this_index_i += 1
                                this_index_j += 1
                            # If we find a bracket in the sequence, we keep the
                            # next three letters and create a new chemical
                            # component which we assign the resepective canonical
                            # code from the canonical sequence
                            elif cur_entity_polymer_one_letter_code[this_index_i] == '(':
                                this_new_chem_comp_name = \
                                    cur_entity_polymer_one_letter_code[this_index_i + 1:this_index_i + 4]
                                this_canonical_res = \
                                    cur_entity_polymer_one_letter_code_canonical[this_index_j]
                                # Create the new chemical component
                                if cur_entity_polymer_type == 'polypeptide(D)':
                                    this_new_chem_comp = \
                                        ihm.DPeptideChemComp(
                                            id=this_new_chem_comp_name,
                                            code=this_new_chem_comp_name,
                                            code_canonical=this_canonical_res)
                                if cur_entity_polymer_type == 'polypeptide(L)':
                                    this_new_chem_comp = \
                                        ihm.LPeptideChemComp(
                                            id=this_new_chem_comp_name,
                                            code=this_new_chem_comp_name,
                                            code_canonical=this_canonical_res)
                                if cur_entity_polymer_type == 'polyribonucleotide':
                                    this_new_chem_comp = \
                                        ihm.RNAChemComp(
                                            id=this_new_chem_comp_name,
                                            code=this_new_chem_comp_name,
                                            code_canonical=this_canonical_res)
                                if cur_entity_polymer_type == 'polydeoxribonucleotide':
                                    this_new_chem_comp = \
                                        ihm.DNAChemComp(
                                            id=this_new_chem_comp_name,
                                            code=this_new_chem_comp_name,
                                            code_canonical=this_canonical_res)
                                list_of_sequence.append(this_new_chem_comp)
                                # And we go to the next entry in the non-canonical
                                # sequence
                                this_index_i += 5
                                this_index_j += 1
                            else:
                                pass
                        # Create the entity with the sequence including the
                        # modified residues.
                        cur_ihm_entity = ihm.Entity(
                            sequence=list_of_sequence,
                            alphabet=cur_alphabet,
                            description=cur_entity_description,
                            source=cur_source)

            # If the entity is a non-polymer
            if cur_entity_type == 'non-polymer':
                cur_entity = ihm.Entity(
                    sequence=[ihm.NonPolymerChemComp(
                        id=cur_entity_nonpolymer_chem_comp_id)],
                    description=cur_ihm_entity_description)
            # If the entity is water
            if cur_entity_type == 'water':
                cur_entity = ihm.Entity(
                    sequence=[ihm.WaterChemComp()],
                    description=cur_ihm_entity_description)
            cur_ihm_entity_index = -1
            if cur_entity not in self.list_ihm_entities:
                self.list_ihm_entities.append(cur_entity)
                self.list_ihm_entity_ids.append(cur_entity_molecular_entity)
            cur_ihm_entity_index = self.list_ihm_entities.index(cur_entity)

            # possibly generate an entity assembly
            cur_entity_assembly_id = \
                xls_data['IHM_Entity_Entity_assembly_id'][i]
            cur_entity_number_of_copies = \
                xls_data['IHM_Entity_Number_of_copies'][i]
            # if the entity_assembly_id is not in the list of entity_assembly_ids
            if cur_entity_assembly_id not in self.list_entity_assembly_ids:
                # create the entity assembly
                cur_entity_assembly = ihm.flr.EntityAssembly()
                # add the id and the entity assembly to the respective lists
                self.list_entity_assembly_ids.append(cur_entity_assembly_id)
                self.list_entity_assemblies.append(cur_entity_assembly)
            # and add the current entity to the respective assembly
            # get the index of the entity_assembly_id
            cur_entity_assembly_id_index = \
                self.list_entity_assembly_ids.index(cur_entity_assembly_id)
            # and add the entity to the entity assembly
            self.list_entity_assemblies[cur_entity_assembly_id_index].add_entity(
                entity=cur_entity,
                num_copies=cur_entity_number_of_copies)

        self.system.entities.extend(self.list_ihm_entities)

    def add_ihm_comparative_starting_model(self,
                                       sheet_name='Starting_comparative_model',
                                       skiprows=3,
                                       header=0):
        """Read and collect information on comparative starting models
        (Templates)"""
        if self.verbose:
            print(" ... Processing tab \'Starting_comparative_model\' ...")
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_starting_comparative_model_Ordinal'])

        for i in range(nr_entries):
            cur_start_model_comp_ordinal = \
                xls_data['IHM_starting_comparative_model_Ordinal'][i]
            cur_start_model_comp_starting_model_id = \
                xls_data['IHM_starting_comparative_model_Starting_model_id'][i]
            cur_start_model_comp_auth_asym_id = \
                xls_data['IHM_starting_comparative_model_Auth_asym_id'][i]
            cur_start_model_comp_seq_id_begin = \
                xls_data['IHM_starting_comparative_model_Seq_id_begin'][i]
            cur_start_model_comp_seq_id_end = \
                xls_data['IHM_starting_comparative_model_Seq_id_end'][i]
            cur_start_model_comp_template_auth_asym_id = \
                xls_data['IHM_starting_comparative_model_Template_auth_asym_id'][i]
            cur_start_model_comp_template_seq_id_begin = \
                xls_data['IHM_starting_comparative_model_Template_seq_id_begin'][i]
            cur_start_model_comp_template_seq_id_end = \
                xls_data['IHM_starting_comparative_model_Template_seq_id_end'][i]
            cur_start_model_comp_template_seq_identity = \
                xls_data['IHM_starting_comparative_model_Template_sequence_identity'][i]
            cur_start_model_comp_template_seq_identity_denominator_id = \
                xls_data['IHM_starting_comparative_model_Template_sequence_identity_denominator_id'][i]
            cur_start_model_comp_dataset_list_id = \
                xls_data['IHM_starting_comparative_model_Dataset_list_id'][i]
            cur_start_model_comp_alignment_file_id = None if (
                    'IHM_starting_comparative_model_Alignment_file_id' not in
                    xls_data.keys() or pandas.isnull(
                        xls_data['IHM_starting_comparative_model_Alignment_file_id'][i]))\
                else xls_data['IHM_starting_comparative_model_Alignment_file_id'][i]

            # translate the sequence identity denominator for the ihm-python implementation
            cur_start_model_comp_template_seq_identity_denominator_id_enum = 0
            if cur_start_model_comp_template_seq_identity_denominator_id == \
                    'Length of the shorter sequence':
                cur_start_model_comp_template_seq_identity_denominator_id_enum = 1
            elif cur_start_model_comp_template_seq_identity_denominator_id == \
                    'Number of aligned positions (including gaps)':
                cur_start_model_comp_template_seq_identity_denominator_id_enum = 2
            elif cur_start_model_comp_template_seq_identity_denominator_id == \
                    'Number of aligned residue pairs (not including gaps)':
                cur_start_model_comp_template_seq_identity_denominator_id_enum = 3
            elif cur_start_model_comp_template_seq_identity_denominator_id == \
                    'Arithmetic mean sequence length':
                cur_start_model_comp_template_seq_identity_denominator_id_enum = 4
            else:  # Other
                cur_start_model_comp_template_seq_identity_denominator_id_enum = 5
            cur_start_model_comp_seq_identity_denominator_object = \
                ihm.startmodel.SequenceIdentityDenominator(
                    cur_start_model_comp_template_seq_identity_denominator_id_enum)
            cur_start_model_comp_seq_identity_object = \
                ihm.startmodel.SequenceIdentity(
                    value=float(
                        cur_start_model_comp_template_seq_identity),
                    denominator=cur_start_model_comp_seq_identity_denominator_object)
            cur_start_model_comp_template_object = \
                ihm.startmodel.Template(
                    dataset=self.list_datasets[self.list_dataset_ids.index(
                        cur_start_model_comp_dataset_list_id)],
                    asym_id=cur_start_model_comp_template_auth_asym_id,
                    seq_id_range=(cur_start_model_comp_seq_id_begin,
                                  cur_start_model_comp_seq_id_end),
                    template_seq_id_range=(cur_start_model_comp_template_seq_id_begin,
                                           cur_start_model_comp_template_seq_id_end),
                    sequence_identity=cur_start_model_comp_seq_identity_object,
                    alignment_file=None if (
                            cur_start_model_comp_alignment_file_id is None)
                    else self.list_external_files_locations[
                        self.list_external_files_ids.index(
                            cur_start_model_comp_alignment_file_id)])

            if not occurs_in_list(cur_start_model_comp_template_object,
                                  self.list_starting_model_templates):
                self.list_starting_model_templates.append(
                    cur_start_model_comp_template_object)
            if cur_start_model_comp_starting_model_id not in \
                    self.tmp_dict_templates_for_starting_models.keys():
                self.tmp_dict_templates_for_starting_models[
                    cur_start_model_comp_starting_model_id] = []
            self.tmp_dict_templates_for_starting_models[
                cur_start_model_comp_starting_model_id].append(
                cur_start_model_comp_template_object)

    def add_ihm_instance_and_AsymUnits(self, sheet_name='Instance', skiprows=3, header=0):
        """Read Instance (AsymUnit) from the excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'Instance\' ...")
        # This part creates the AsymUnit (or AsymUnitRange),
        # ihm_model_representation, and ihm_starting_model_details
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_Instance_Ordinal'])

        for i in range(nr_entries):
            cur_instance_ordinal = xls_data['IHM_Instance_Ordinal']
            cur_instance_chain_id = xls_data['IHM_Instance_Chain_id'][i]
            cur_instance_entity_id = xls_data['IHM_Instance_Entity_id'][i]
            cur_instance_details = \
                None if ('IHM_Instance_Details' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Instance_Details'][i])) \
                    else xls_data['IHM_Instance_Details'][i]
            cur_instance_seq_begin = \
                None if ('IHM_Instance_Seq_begin' not in xls_data.keys()
                         or pandas.isnull(
                            xls_data['IHM_Instance_Seq_begin'][i])) \
                    else int(xls_data['IHM_Instance_Seq_begin'][i])
            cur_instance_seq_end = \
                None if ('IHM_Instance_Seq_end' not in xls_data.keys() or
                         pandas.isnull(xls_data['IHM_Instance_Seq_end'][i])) \
                    else int(xls_data['IHM_Instance_Seq_end'][i])
            cur_instance_structure_assembly_id = \
                xls_data['IHM_Instance_Structure_assembly_id'][i]
            cur_instance_structure_assembly_name = \
                xls_data['IHM_Instance_Structure_assembly_name'][i]
            cur_instance_structure_assembly_description = \
                xls_data['IHM_Instance_Structure_assembly_description'][i]
            cur_instance_model_representation_id = \
                None if ('IHM_Instance_Model_representation_id' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Model_representation_id'][i])) \
                    else xls_data['IHM_Instance_Model_representation_id'][i]
            cur_instance_model_object_primitive = \
                None if ('IHM_Instance_Model_object_primitive' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Model_object_primitive'][i])) \
                    else xls_data['IHM_Instance_Model_object_primitive'][i]
            # True corresponds to rigid, while False corresponds to flexible
            cur_instance_model_mode = \
                None if ('IHM_Instance_Model_mode' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Model_mode'][i])) \
                    else (True if
                          (xls_data['IHM_Instance_Model_mode'][i] == 'rigid')
                          else False)
            cur_instance_model_granularity = \
                None if ('IHM_Instance_Model_granularity' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Model_granularity'][i])) \
                    else xls_data['IHM_Instance_Model_granularity'][i]
            cur_instance_object_count = \
                None if ('IHM_Instance_Object_count' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Object_count'][i])) \
                    else xls_data['IHM_Instance_Object_count'][i]
            cur_instance_starting_model_id = \
                None if ('IHM_Instance_Starting_model_id' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Starting_model_id'][i])) \
                    else xls_data['IHM_Instance_Starting_model_id'][i]
            cur_instance_starting_model_source = \
                None if ('IHM_Instance_Starting_model_source' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Starting_model_source'][i])) \
                    else xls_data['IHM_Instance_Starting_model_source'][i]
            cur_instance_starting_model_chain_id = \
                xls_data['IHM_Instance_Starting_model_chain_id'][i]
            cur_instance_starting_model_sequence_offset = \
                None if ('IHM_Instance_Starting_model_sequence_offset' not in
                         xls_data.keys()) else \
                    (0 if (pandas.isnull(
                        xls_data['IHM_Instance_Starting_model_sequence_offset'][i]))
                     else int(xls_data['IHM_Instance_Starting_model_sequence_offset'][i]))
            cur_instance_starting_model_dataset_list_id = \
                xls_data['IHM_Instance_Starting_model_dataset_list_id'][i]
            # Use external_file_id if available
            cur_instance_starting_model_external_file_id = \
                None if ('IHM_Instance_Starting_model_external_file_id' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Instance_Starting_model_external_file_id'][i])) \
                    else xls_data['IHM_Instance_Starting_model_external_file_id'][i]

            # Asym Units
            # get the respective entity
            cur_entity_for_index = \
                self.list_ihm_entities[
                    self.list_ihm_entity_ids.index(cur_instance_entity_id)]
            cur_asym_unit = ihm.AsymUnit(
                entity=cur_entity_for_index,
                details=cur_instance_details,
                id=cur_instance_chain_id)
            # check whether there is already an asymmetric unit with the same id
            if cur_instance_chain_id not in self.list_asym_units_ids:
                # Comparing the asym_units directly does not work, because the
                # details are different. This would result in later problems
                # when  the dumper finds the same asym_unit_id twice
#             if (not occurs_in_list(cur_asym_unit, list_asym_units)) and
#             (cur_ihm_instance_chain_id not in list_asym_units_ids):
                self.list_asym_units.append(cur_asym_unit)
                self.list_asym_units_ids.append(cur_instance_chain_id)
            else:
                # otherwise, we use the asym_unit that was created before
#            cur_asym_unit = [x for x in list_asym_units if x.__dict__ == cur_asym_unit.__dict__][0]
                cur_asym_unit = \
                    self.list_asym_units[
                        self.list_asym_units_ids.index(cur_instance_chain_id)]
                ## join the details for the Asym-units
                cur_asym_unit.details = \
                    cur_asym_unit.details if cur_instance_details in \
                                             cur_asym_unit.details \
                        else " & ".join([cur_asym_unit.details,
                                         cur_instance_details])

            # Asym Unit Ranges
            cur_asym_unit_range = ihm.AsymUnitRange(
                asym=self.list_asym_units[self.list_asym_units.index(cur_asym_unit)],
                seq_id_begin=cur_instance_seq_begin,
                seq_id_end=cur_instance_seq_end)
            if not occurs_in_list(cur_asym_unit_range,
                                  self.list_asym_unit_ranges):
                self.list_asym_unit_ranges.append(cur_asym_unit_range)

            # Starting model details
            cur_starting_model = None
            cur_templates = None
            if cur_instance_starting_model_id is not None:
                if cur_instance_starting_model_source == 'comparative model':
                    cur_templates = \
                        self.tmp_list_templates_for_starting_models[
                            cur_instance_starting_model_id]
                # use the external file dataset list id if available
                if cur_instance_starting_model_external_file_id is not None:
                    cur_starting_model = ihm.startmodel.StartingModel(
                        asym_unit=self.list_asym_unit_ranges[
                            self.list_asym_unit_ranges.index(cur_asym_unit_range)],
                        dataset=self.list_connection_dataset_external_file_datasets[
                            list_connection_dataset_external_file_ids.index(
                                cur_instance_starting_model_external_file_id)],
                        asym_id=cur_instance_starting_model_chain_id,
                        offset=cur_instance_starting_model_sequence_offset,
                        templates=cur_templates)
                else:
                    cur_starting_model = ihm.startmodel.StartingModel(
                        self.list_asym_unit_ranges[
                            self.list_asym_unit_ranges.index(cur_asym_unit_range)],
                        dataset=self.list_datasets[
                            self.list_dataset_ids.index(
                                cur_instance_starting_model_dataset_list_id)],
                        asym_id=cur_instance_starting_model_chain_id,
                        offset=cur_instance_starting_model_sequence_offset,
                        templates=cur_templates)

            # Model representation
            # !!! Note: Only by-atom tested so far
            if cur_instance_model_representation_id is not None:
                if cur_instance_model_granularity == 'by-atom':
                    cur_model_representation = \
                        ihm.representation.AtomicSegment(
                            asym_unit=cur_asym_unit_range,
                            rigid=cur_instance_model_mode,
                            starting_model=cur_starting_model)
                if cur_instance_model_granularity == 'by-residue':
                    cur_model_representation = \
                        ihm.representation.ResidueSegment(
                            asym_unit=cur_asym_unit_range,
                            rigid=cur_instance_model_mode,
                            primitive=cur_instance_model_object_primitive,
                            starting_model=cur_starting_model)
                if cur_instance_model_granularity == 'multi-residue':
                    cur_model_representation = \
                        ihm.representation.MultiResidueSegment(
                            asym_unit=cur_asym_unit_range,
                            rigid=cur_instance_model_mode,
                            primitive=cur_instance_model_object_primitive,
                            starting_model=cur_starting_model)
                if cur_instance_model_granularity == 'by-feature':
                    cur_model_representation = \
                        ihm.representation.FeatureSegment(
                            asym_unit=cur_asym_unit_range,
                            rigid=cur_instance_model_mode,
                            primitive=cur_instance_model_object_primitive,
                            count=cur_instance_object_count,
                            starting_model=cur_starting_model)
                if cur_instance_model_representation_id not in \
                        self.tmp_dict_for_model_representations.keys():
                    self.tmp_dict_for_model_representations[
                        cur_instance_model_representation_id] = []
                self.tmp_dict_for_model_representations[
                    cur_instance_model_representation_id].append(
                    cur_model_representation)

            # structure assembly by id
            if cur_instance_structure_assembly_id not in \
                    self.tmp_dict_for_structure_assembly.keys():
                self.tmp_dict_for_structure_assembly[
                    cur_instance_structure_assembly_id] = []
                self.tmp_dict_for_structure_assembly_names[
                    cur_instance_structure_assembly_id] = \
                    (cur_instance_structure_assembly_name)
                self.tmp_dict_for_structure_assembly_description[
                    cur_instance_structure_assembly_id] = \
                    (cur_instance_structure_assembly_description)

            if cur_asym_unit not in self.tmp_dict_for_structure_assembly[
                cur_instance_structure_assembly_id]:
                self.tmp_dict_for_structure_assembly[
                    cur_instance_structure_assembly_id].append(cur_asym_unit)
            self.tmp_dict_for_structure_assembly_names[
                cur_instance_structure_assembly_id] = (
                cur_instance_structure_assembly_name)
            self.tmp_dict_for_structure_assembly_description[
                cur_instance_structure_assembly_id] = (
                cur_instance_structure_assembly_description)

        # create the model representations
        ### TODO: CHECK: Should this be self.tmp_dict_for_model_representations.keys() ???
        for groupkey in self.tmp_dict_for_model_representations:
            cur_model_representation = \
                ihm.representation.Representation(
                    self.tmp_dict_for_model_representations[groupkey])
            self.list_ihm_model_representations.append(cur_model_representation)
            self.list_ihm_model_representation_ids.append(groupkey)

        ## create the structure assemblies
        for groupkey in self.tmp_dict_for_structure_assembly.keys():
            cur_structure_assembly = ihm.Assembly(
                elements=self.tmp_dict_for_structure_assembly[groupkey],
                name=self.tmp_dict_for_structure_assembly_names[groupkey],
                description=self.tmp_dict_for_structure_assembly_description[groupkey])
            self.list_structure_assemblies.append(cur_structure_assembly)
            self.list_structure_assembly_ids.append(groupkey)

        self.system.asym_units.extend(self.list_asym_units)

    def add_ihm_modeling_protocol(self,
                                  sheet_name_protocol='Modeling_protocol',
                                  skiprows_protocol=3,
                                  header_protocol=0,
                                  sheet_name_post_process='Modeling_post_process',
                                  skiprows_post_process=3,
                                  header_post_process=0):
        """Read the modeling protocol and modeling post process information
        from the excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'Modeling_post_process\' ...")
        xls_data_post = pandas.read_excel(self.xls_file,
                                          sheet_name=sheet_name_post_process,
                                          skiprows=skiprows_post_process,
                                          header=header_post_process)
        nr_entries_post = len(xls_data_post['IHM_Post_process_ID'])

        # Store the modeling analysis steps for the creation of the Analysis object
        tmp_dict_for_protocols_analysis = {}
        tmp_dict_for_protocols_analysis_ids = {}

        for i in range(nr_entries_post):
            cur_post_process_id = \
                xls_data_post['IHM_Post_process_ID'][i]
            cur_post_process_modeling_protocol_id = \
                xls_data_post['IHM_Post_process_Protocol_id'][i]
            cur_post_process_analysis_id = \
                xls_data_post['IHM_Post_process_Step_id'][i]
            cur_post_process_struct_assembly_id = \
                xls_data_post['IHM_Post_process_Structure_assembly_id'][i]
            cur_post_process_dataset_group_id = \
                xls_data_post['IHM_Post_process_Dataset_group_id'][i]
            cur_post_process_type = \
                xls_data_post['IHM_Post_process_Type'][i]
            cur_post_process_feature = \
                xls_data_post['IHM_Post_process_Feature'][i]
            cur_post_process_feature_name = \
                xls_data_post['IHM_Post_process_Feature_name'][i]
            cur_post_process_number_model_begin = \
                None if ('IHM_Post_process_Number_model_begin' not in
                         xls_data_post.keys() or pandas.isnull(
                            xls_data_post['IHM_Post_process_Number_model_begin'][i])) \
                    else int(xls_data_post['IHM_Post_process_Number_model_begin'][i])
            cur_post_process_number_model_end = \
                None if ('IHM_Post_process_Number_model_end' not in
                         xls_data_post.keys() or pandas.isnull(
                            xls_data_post['IHM_Post_process_Number_model_end'][i])) \
                    else int(xls_data_post['IHM_Post_process_Number_model_end'][i])
            cur_post_process_software_id = \
                xls_data_post['IHM_Post_process_Software_id'][i]

            # Post-processing steps (ihm-analysis)
            cur_step = None
            if cur_post_process_type == 'cluster':
                cur_step = ihm.analysis.ClusterStep(
                    feature=cur_post_process_feature,
                    num_models_begin=cur_post_process_number_model_begin,
                    num_models_end=cur_post_process_number_model_end,
                    assembly=self.list_structure_assemblies[
                        self.list_structure_assembly_ids.index(
                            cur_post_process_struct_assembly_id)],
                    dataset_group=self.list_dataset_groups[
                        self.list_dataset_group_ids.index(
                            cur_post_process_dataset_group_id)],
                    software=self.list_ihm_softwares[
                        self.list_ihm_software_ids.index(
                            cur_post_process_software_id)])
            if cur_post_process_type == 'filter':
                cur_step = ihm.analysis.FilterStep(
                    feature=cur_post_process_feature,
                    num_models_begin=cur_post_process_number_model_begin,
                    num_models_end=cur_post_process_number_model_end,
                    assembly=self.list_structure_assemblies[
                        self.list_structure_assembly_ids.index(
                            cur_post_process_struct_assembly_id)],
                    dataset_group=self.list_dataset_groups[
                        self.list_dataset_group_ids.index(
                            cur_post_process_dataset_group_id)],
                    software=self.list_ihm_softwares[
                        self.list_ihm_software_ids.index(
                            cur_post_process_software_id)])

            if cur_post_process_type == 'rescore':
                cur_step = ihm.analysis.RescoreStep(
                    feature=cur_post_process_feature,
                    num_models_begin=cur_post_process_number_model_begin,
                    num_models_end=cur_post_process_number_model_end,
                    assembly=self.list_structure_assemblies[
                        self.list_structure_assembly_ids.index(
                            cur_post_process_struct_assembly_id)],
                    dataset_group=self.list_dataset_groups[
                        self.list_dataset_group_ids.index(
                            cur_post_process_dataset_group_id)],
                    software=self.list_ihm_softwares[
                        self.list_ihm_software_ids.index(
                            cur_post_process_software_id)])

            if cur_post_process_type == 'validation':
                cur_step = ihm.analysis.ValidationStep(
                    feature=cur_post_process_feature,
                    num_models_begin=cur_post_process_number_model_begin,
                    num_models_end=cur_post_process_number_model_end,
                    assembly=self.list_structure_assemblies[
                        self.list_structure_assembly_ids.index(
                            cur_post_process_struct_assembly_id)],
                    dataset_group=self.list_dataset_groups[
                        self.list_dataset_group_ids.index(
                            cur_post_process_dataset_group_id)],
                    software=self.list_ihm_softwares[
                        self.list_ihm_software_ids.index(
                            cur_post_process_software_id)])

            if cur_post_process_type == 'none':
                cur_step = ihm.analysis.EmptyStep(
                    feature=cur_post_process_feature,
                    num_models_begin=cur_post_process_number_model_begin,
                    num_models_end=cur_post_process_number_model_end,
                    assembly=self.list_structure_assemblies[
                        self.list_structure_assembly_ids.index(
                            cur_post_process_struct_assembly_id)],
                    dataset_group=self.list_dataset_groups[
                        self.list_dataset_group_ids.index(
                            cur_post_process_dataset_group_id)],
                    software=self.list_ihm_softwares[
                        self.list_ihm_software_ids.index(
                            cur_post_process_software_id)])

            if cur_post_process_type == 'other':
                cur_step = ihm.analysis.Step(
                    feature=cur_post_process_feature,
                    num_models_begin=cur_post_process_number_model_begin,
                    num_models_end=cur_post_process_number_model_end,
                    assembly=self.list_structure_assemblies[
                        self.list_structure_assembly_ids.index(
                            cur_post_process_struct_assembly_id)],
                    dataset_group=self.list_dataset_groups[
                        self.list_dataset_group_ids.index(
                            cur_post_process_dataset_group_id)],
                    software=self.list_ihm_softwares[
                        self.list_ihm_software_ids.index(
                            cur_post_process_software_id)])

            # store the step for later addition to an Analysis object and also
            # store the id
            if cur_step not in self.list_ihm_modeling_protocol_analysis_steps:
                self.list_ihm_modeling_protocol_analysis_steps.append(cur_step)
                self.list_ihm_modeling_protocol_analysis_step_ids.append(
                    cur_post_process_id)
            # possibly create a new list for one Analysis to save the steps to
            if cur_post_process_modeling_protocol_id not in \
                    tmp_dict_for_protocols_analysis.keys():
                tmp_dict_for_protocols_analysis[
                    cur_post_process_modeling_protocol_id] = []
                tmp_dict_for_protocols_analysis_ids[
                    cur_post_process_modeling_protocol_id] = []
            # and save the steps to the respective analysis id
            tmp_dict_for_protocols_analysis[
                cur_post_process_modeling_protocol_id].append(cur_step)
            tmp_dict_for_protocols_analysis_ids[
                cur_post_process_modeling_protocol_id].append(
                cur_post_process_id)

        if self.verbose:
            print(' ... Processing tab \'Modeling_protocol\' ...')
        xls_data_protocol = pandas.read_excel(self.xls_file,
                                              sheet_name='Modeling_protocol',
                                              skiprows=3,
                                              header=0)
        nr_entries_protocol = \
            len(xls_data_protocol['IHM_Modeling_protocol_Ordinal'])

        ### !! TODO: check whether these dicts are only needed in this function
        tmp_dict_for_protocols_modeling = {}
        tmp_dict_for_protocols_modeling_ids = {}
        tmp_dict_for_protocols_names = {}
        tmp_dict_store_analyses_for_protocols = {}

        for i in range(nr_entries_protocol):
            cur_modeling_protocol_ordinal = \
                xls_data_protocol['IHM_Modeling_protocol_Ordinal'][i]
            cur_modeling_protocol_protocol_id = \
                xls_data_protocol['IHM_Modeling_protocol_Protocol_id'][i]
            cur_modeling_protocol_step_id = \
                xls_data_protocol['IHM_Modeling_protocol_Step_id'][i]
            cur_modeling_protocol_structure_assembly_id = \
                xls_data_protocol['IHM_Modeling_protocol_Structure_assembly_id'][i]
            cur_modeling_protocol_dataset_group_id = \
                xls_data_protocol['IHM_Modeling_protocol_Dataset_group_id'][i]
            cur_modeling_protocol_protocol_name = \
                xls_data_protocol['IHM_Modeling_protocol_Protocol_name'][i]
            cur_modeling_protocol_step_name = \
                xls_data_protocol['IHM_Modeling_protocol_Step_name'][i]
            cur_modeling_protocol_step_method = \
                None if ('IHM_Modeling_protocol_Step_method' not in
                         xls_data_protocol.keys() or pandas.isnull(
                            xls_data_protocol['IHM_Modeling_protocol_Step_method'][i])) \
                    else xls_data_protocol['IHM_Modeling_protocol_Step_method'][i]
            cur_modeling_protocol_number_model_begin = \
                None if ('IHM_Modeling_protocol_Number_model_begin' not in
                         xls_data_protocol.keys() or pandas.isnull(
                            xls_data_protocol['IHM_Modeling_protocol_Number_model_begin'][i])) \
                    else int(xls_data_protocol['IHM_Modeling_protocol_Number_model_begin'][i])
            cur_modeling_protocol_number_model_end = \
                None if ('IHM_Modeling_protocol_Number_model_end' not in
                         xls_data_protocol.keys() or pandas.isnull(
                            xls_data_protocol['IHM_Modeling_protocol_Number_model_end'][i])) \
                    else int(xls_data_protocol['IHM_Modeling_protocol_Number_model_end'][i])
            cur_modeling_protocol_multi_scale_flag = \
                xls_data_protocol['IHM_Modeling_protocol_Multi_scale_flag'][i] in \
                ['Yes', 'YES', 'yes']
            cur_modeling_protocol_multi_state_flag = \
                xls_data_protocol['IHM_Modeling_protocol_Multi_state_flag'][i] in \
                ['Yes', 'YES', 'yes']
            cur_modeling_protocol_ordered_flag = \
                xls_data_protocol['IHM_Modeling_protocol_Ordered_flag'][i] in \
                ['Yes', 'YES', 'yes']
            cur_modeling_protocol_ensemble_flag = \
                xls_data_protocol['IHM_Modeling_protocol_Ensemble_flag'][i] in \
                ['Yes', 'YES', 'yes']
            cur_modeling_protocol_software_id = \
                xls_data_protocol['IHM_Modeling_protocol_Software_id'][i]

            # store the software used for the current modeling protocol ordinal
            self.dict_ihm_modeling_protocol_software_ids[
                cur_modeling_protocol_ordinal] = \
                cur_modeling_protocol_software_id

            cur_step = None
            # Protocol steps
            cur_step = ihm.protocol.Step(
                assembly=self.list_structure_assemblies[
                    self.list_structure_assembly_ids.index(
                        cur_modeling_protocol_structure_assembly_id)],
                dataset_group=self.list_dataset_groups[
                    self.list_dataset_group_ids.index(
                        cur_modeling_protocol_dataset_group_id)],
                method=cur_modeling_protocol_step_method,
                num_models_begin=cur_modeling_protocol_number_model_begin,
                num_models_end=cur_modeling_protocol_number_model_end,
                software=self.list_ihm_softwares[
                    self.list_ihm_software_ids.index(
                        cur_modeling_protocol_software_id)],
                multi_scale=cur_modeling_protocol_multi_scale_flag,
                multi_state=cur_modeling_protocol_multi_state_flag,
                ordered=cur_modeling_protocol_ordered_flag,
                name=cur_modeling_protocol_step_name)
            # save the current step for the modeling protocols
            if cur_step not in self.list_ihm_modeling_protocol_modeling_steps:
                self.list_ihm_modeling_protocol_modeling_steps.append(cur_step)
                self.list_ihm_modeling_protocol_modeling_step_ids.append(
                    cur_modeling_protocol_step_id)
            # create the entry for the current protocol_id
            if cur_modeling_protocol_protocol_id not in \
                    tmp_dict_for_protocols_modeling.keys():
                tmp_dict_for_protocols_modeling[
                    cur_modeling_protocol_protocol_id] = []
                tmp_dict_for_protocols_modeling_ids[
                    cur_modeling_protocol_protocol_id] = []
            # and store the entries to the current protocol id
            tmp_dict_for_protocols_modeling[
                cur_modeling_protocol_protocol_id].append(cur_step)
            tmp_dict_for_protocols_modeling_ids[
                cur_modeling_protocol_protocol_id].append(
                cur_modeling_protocol_step_id)
            tmp_dict_for_protocols_names[
                cur_modeling_protocol_protocol_id] = \
                cur_modeling_protocol_protocol_name

        # Analyses performed
        for groupkey in tmp_dict_for_protocols_analysis.keys():
            # create an Analysis object
            cur_modeling_protocol_analysis = ihm.analysis.Analysis()
            # store each analysis step that was assigned to the current analysis
            # id (by its cur_ihm_modeling_protocol_post_process_analysis_id)
            for entry in tmp_dict_for_protocols_analysis[groupkey]:
                cur_modeling_protocol_analysis.steps.append(entry)
            if cur_modeling_protocol_analysis not in \
                    self.list_ihm_modeling_protocol_analyses:
                self.list_ihm_modeling_protocol_analyses.append(
                    cur_modeling_protocol_analysis)
                self.list_ihm_modeling_protocol_analyses_ids.append(groupkey)

        # Protocol
        for groupkey in tmp_dict_for_protocols_modeling.keys():
            # create a protocol with the respective key
            cur_modeling_protocol = ihm.protocol.Protocol(
                name=tmp_dict_for_protocols_names[groupkey])
            for entry in tmp_dict_for_protocols_modeling[groupkey]:
                cur_modeling_protocol.steps.append(entry)

            # get the analyses_id with the respective groupkey
            # only if there were post-processing steps
            if len(self.list_ihm_modeling_protocol_analyses_ids) > 0:
                cur_analysis_id_index = \
                    self.list_ihm_modeling_protocol_analyses_ids.index(groupkey)
                cur_modeling_protocol.analyses.append(
                    self.list_ihm_modeling_protocol_analyses[
                        cur_analysis_id_index])

            # finally, add the modeling protocol to the list of modeling protocols
            if cur_modeling_protocol not in \
                    self.list_ihm_modeling_protocols:
                self.list_ihm_modeling_protocols.append(cur_modeling_protocol)
                self.list_ihm_modeling_protocols_ids.append(groupkey)



    def add_ihm_multi_state_modeling_and_models(
            self,
            sheet_name_multi_state_modeling='Multi-state modeling',
            skiprows_multi_state_modeling=3,
            header_multi_state_modeling=0,
            sheet_name_models='Models',
            skiprows_models=3,
            header_models=0):
        """Read multi-state modeling information from the excel sheet"""
        # Multi-state modeling
        if self.verbose:
            print(" ... Processing tab \'Multi-state modeling\' ...")
        xls_data_multi = pandas.read_excel(
            self.xls_file,
            sheet_name=sheet_name_multi_state_modeling,
            skiprows=skiprows_multi_state_modeling,
            header=header_multi_state_modeling)
        nr_entries_multi = \
            len(xls_data_multi['IHM_Multi_state_modeling_Ordinal_id'])

        tmp_dict_collect_for_multi_state = {}
        for i in range(nr_entries_multi):
            cur_multi_state_modeling_ordinal_id = \
                xls_data_multi['IHM_Multi_state_modeling_Ordinal_id'][i]
            cur_multi_state_modeling_state_id = \
                xls_data_multi['IHM_Multi_state_modeling_state_id'][i]
            cur_multi_state_modeling_state_group_id = \
                xls_data_multi['IHM_Multi_state_modeling_state_group_id'][i]
            cur_multi_state_modeling_population_fraction = \
                xls_data_multi['IHM_Multi_state_modeling_population_fraction'][i]
            cur_multi_state_modeling_population_fraction_sd = \
                xls_data_multi['IHM_Multi_state_modeling_population_fraction_standard_deviation'][i]
            cur_multi_state_modeling_state_type = \
                xls_data_multi['IHM_Multi_state_modeling_state_type'][i]
            cur_multi_state_modeling_state_name = \
                xls_data_multi['IHM_Multi_state_modeling_state_name'][i]
            cur_multi_state_modeling_state_details = \
                None if ('IHM_Multi_state_modeling_state_details' not in
                         xls_data_multi.keys() or pandas.isnull(
                            xls_data_multi['IHM_Multi_state_modeling_state_details'][i])) \
                    else xls_data_multi['IHM_Multi_state_modeling_state_details'][i]
            cur_multi_state_modeling_model_group_id = \
                None if ('IHM_Multi_state_modeling_model_group_id' not in
                         xls_data_multi.keys() or pandas.isnull(
                            xls_data_multi['IHM_Multi_state_modeling_model_group_id'][i])) \
                    else xls_data_multi['IHM_Multi_state_modeling_model_group_id'][i]
            cur_multi_state_modeling_experiment_type = \
                xls_data_multi['IHM_Multi_state_modeling_experiment_type'][i]

            # If the current state_id is not an entry in tmp_list_collect_for_multi_state yet
            if cur_multi_state_modeling_state_id not in \
                    tmp_dict_collect_for_multi_state.keys():
                tmp_dict_collect_for_multi_state[
                    cur_multi_state_modeling_state_id] = {}
            tmp_dict_collect_for_multi_state[
                cur_multi_state_modeling_state_id]['state_type'] = \
                cur_multi_state_modeling_state_type
            tmp_dict_collect_for_multi_state[
                cur_multi_state_modeling_state_id]['state_name'] = \
                cur_multi_state_modeling_state_name
            tmp_dict_collect_for_multi_state[
                cur_multi_state_modeling_state_id]['state_details'] = \
                cur_multi_state_modeling_state_details
            tmp_dict_collect_for_multi_state[
                cur_multi_state_modeling_state_id]['experiment_type'] = \
                cur_multi_state_modeling_experiment_type
            tmp_dict_collect_for_multi_state[
                cur_multi_state_modeling_state_id]['population_fraction'] = \
                cur_multi_state_modeling_population_fraction

            # If there is no model_group_ID, then the state will be added without
            # a model group
            # This is to allow states, for which distance restraints are derived,
            # but no structural models are deposited
            if cur_multi_state_modeling_model_group_id is None:
                #            dummy_model_group = ihm.model.ModelGroup(elements = [])
                cur_state = ihm.model.State(
                    elements=[],
                    type=cur_multi_state_modeling_state_type,
                    name=cur_multi_state_modeling_state_name,
                    details=cur_multi_state_modeling_state_details,
                    experiment_type=cur_multi_state_modeling_experiment_type,
                    population_fraction=cur_multi_state_modeling_population_fraction)

                if cur_state not in list_empty_states:
                    self.list_empty_states.append(cur_state)
                    self.list_empty_state_ids.append(cur_multi_state_modeling_state_id)

        # Models
        if self.verbose:
            print(' ... Processing tab \'Models\' ...')
        xls_ihm_models = pandas.read_excel(
            self.xls_file,
            sheet_name=sheet_name_models,
            skiprows=skiprows_models,
            header=header_models)
        nr_entries_models = len(xls_ihm_models['IHM_Models_Model_number'])

        tmp_dict_for_model_group_models = {}
        tmp_dict_for_model_group_names = {}
        tmp_dict_for_multi_state_models = {}
        tmp_dict_for_multi_state_names = {}
        # collect model groups for the connection to the multi_state_models
        tmp_dict_multi_state_modeling_model_group_link = {}

        for i in range(nr_entries_models):
            cur_models_model_number = \
                xls_ihm_models['IHM_Models_Model_number'][i]
            cur_models_model_name = \
                xls_ihm_models['IHM_Models_Model_name'][i]
            cur_models_model_group_id = \
                xls_ihm_models['IHM_Models_Model_group_id'][i]
            cur_models_model_group_name = \
                xls_ihm_models['IHM_Models_Model_group_name'][i]
            cur_models_multi_state_id = \
                xls_ihm_models['IHM_Models_Multi_state_id'][i]
            # either get the value from here, or get it automatically from the
            # dictionary tmp_list_collect_for_multi_state by the
            # cur_ihm_models_multi_state_id
            # tmp_list_collect_for_multi_state[cur_ihm_models_multi_state_id]['state_name']
            cur_models_multi_state_name = \
                xls_ihm_models['IHM_Models_Multi_state_name'][i]
            cur_models_assembly_id = \
                xls_ihm_models['IHM_Models_assembly_id'][i]
            cur_models_representation_id = \
                xls_ihm_models['IHM_Models_Representation_id'][i]
            cur_models_protocol_id = \
                xls_ihm_models['IHM_Models_Protocol_id'][i]
            cur_models_model_representative = \
                xls_ihm_models['IHM_Models_Model_representative'][i]
            cur_model_representative_selection_criteria = \
                None if ('IHM_Models_Model_representative_selection_criteria' not in
                         xls_ihm_models.keys() or pandas.isnull(
                            xls_ihm_models['IHM_Models_Model_representative_selection_criteria'][i]))\
                    else xls_ihm_models['IHM_Models_Model_representative_selection_criteria'][i]

            cur_model = ihm.model.Model(
                assembly=self.list_structure_assemblies[
                    self.list_structure_assembly_ids.index(
                        cur_models_assembly_id)],
                protocol=self.list_ihm_modeling_protocols[
                    self.list_ihm_modeling_protocols_ids.index(
                        cur_models_protocol_id)],
                representation=self.list_ihm_model_representations[
                    self.list_ihm_model_representation_ids.index(
                        cur_models_representation_id)],
                name=cur_models_model_name)
            if cur_model not in self.list_models:
                self.list_models.append(cur_model)
                self.list_models_ids.append(cur_models_model_number)

            # store the information for the model groups
            # (those will be created once all the information is collected)
            if cur_models_model_group_id not in \
                    tmp_dict_for_model_group_models.keys():
                tmp_dict_for_model_group_models[cur_models_model_group_id] = []
            tmp_dict_for_model_group_models[
                cur_models_model_group_id].append(cur_model)
            tmp_dict_for_model_group_names[cur_models_model_group_id] = \
                cur_models_model_group_name
            if cur_models_multi_state_id not in \
                    tmp_dict_for_multi_state_models.keys():
                tmp_dict_for_multi_state_models[cur_models_multi_state_id] = []
            tmp_dict_for_multi_state_models[
                cur_models_multi_state_id].append(cur_model)
            tmp_dict_for_multi_state_names[cur_models_multi_state_id] = \
                cur_models_multi_state_name
            # store the ID of the model groups for the multi-state-modeling
            if cur_models_multi_state_id not in \
                    tmp_dict_multi_state_modeling_model_group_link:
                tmp_dict_multi_state_modeling_model_group_link[
                    cur_models_multi_state_id] = []
            if cur_models_model_group_id not in \
                    tmp_dict_multi_state_modeling_model_group_link[
                        cur_models_multi_state_id]:
                tmp_dict_multi_state_modeling_model_group_link[
                    cur_models_multi_state_id].append(cur_models_model_group_id)

        # create the model groups (model.py)
        for groupkey in tmp_dict_for_model_group_models.keys():
            cur_model_group = ihm.model.ModelGroup(
                elements=tmp_dict_for_model_group_models[groupkey],
                name=tmp_dict_for_model_group_names[groupkey])
            if cur_model_group not in self.list_model_groups:
                self.list_model_groups.append(cur_model_group)
                self.list_model_group_ids.append(groupkey)

        # create the state (model.py)
        for groupkey in tmp_dict_multi_state_modeling_model_group_link.keys():
            # collect all model groups for the groupkey
            cur_list_of_model_groups = []
            for entry in tmp_dict_multi_state_modeling_model_group_link[groupkey]:
                cur_list_of_model_groups.append(
                    self.list_model_groups[
                        self.list_model_group_ids.index(entry)])
            cur_state = ihm.model.State(
                elements=cur_list_of_model_groups,
                name=tmp_dict_for_multi_state_names[groupkey],
                experiment_type=tmp_dict_collect_for_multi_state[groupkey][
                    'experiment_type'],
                population_fraction=tmp_dict_collect_for_multi_state[groupkey][
                    'population_fraction'],
                type=tmp_dict_collect_for_multi_state[groupkey]['state_type'],
                details=tmp_dict_collect_for_multi_state[groupkey]['state_details'])
            if cur_state not in self.list_models_states:
                self.list_models_states.append(cur_state)
                self.list_models_state_ids.append(groupkey)

        ## TODO? How to use this properly?
        ## create the state group (model.py)
        ## collection of all models
        ## combine all states that include models and those that don't
        self.list_models_states.extend(self.list_empty_states)
        self.list_models_state_ids.extend(self.list_empty_state_ids)
        cur_state_group = ihm.model.StateGroup(elements=self.list_models_states)
        self.system.state_groups.append(cur_state_group)

    def add_ihm_ensemble(self,
                         sheet_name='Ensemble Info',
                         skiprows=3,
                         header=0):
        """Read the IHM ensemble information from the excel sheet."""
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['IHM_Ensemble_ID'])

        for i in range(nr_entries):
            cur_ensemble_id = xls_data['IHM_Ensemble_ID'][i]
            cur_ensemble_name = xls_data['IHM_Ensemble_Name'][i]
            cur_ensemble_post_process_id = \
                xls_data['IHM_Ensemble_Post_process_id'][i]
            cur_ensemble_model_group_id = \
                xls_data['IHM_Ensemble_Model_group_id'][i]
            cur_ensemble_clustering_method = \
                xls_data['IHM_Ensemble_Clustering_method'][i]
            cur_ensemble_clustering_feature = \
                xls_data['IHM_Ensemble_Clustering_feature'][i]
            cur_ensemble_number_of_models_in_ensemble = \
                xls_data['IHM_Ensemble_Number_of_models_in_ensemble'][i]
            cur_ensemble_number_of_models_deposited = \
                xls_data['IHM_Ensemble_Number_of_models_deposited'][i]
            cur_ensemble_precision_value = \
                None if ('IHM_Ensemble_Precision_value' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['IHM_Ensemble_Precision_value'][i])) \
                    else xls_data['IHM_Ensemble_Precision_value'][i]
            cur_ensemble_file_id = \
                None if ('IHM_Ensemble_File_id' not in xls_data.keys() or
                         pandas.isnull(
                             xls_data['IHM_Ensemble_File_id'][i])) \
                    else xls_data['IHM_Ensemble_File_id'][i]

            cur_ihm_ensemble = ihm.model.Ensemble(
                model_group=self.list_model_groups[
                    self.list_model_group_ids.index(
                        cur_ensemble_model_group_id)],
                num_models=cur_ensemble_number_of_models_in_ensemble,
                post_process=self.list_ihm_modeling_protocol_analysis_steps[
                    self.list_ihm_modeling_protocol_analysis_step_ids.index(
                        cur_ensemble_post_process_id)],
                clustering_method=cur_ensemble_clustering_method,
                clustering_feature=cur_ensemble_clustering_feature,
                name=cur_ensemble_name,
                precision=cur_ensemble_precision_value,
                file=None if (cur_ensemble_file_id is None) else
                self.list_external_files_locations[
                    self.list_external_files_ids.index(cur_ensemble_file_id)])

            if not occurs_in_list(cur_ihm_ensemble, self.list_ensembles):
                self.list_ensembles.append(cur_ihm_ensemble)

        # add the ensembles to the system
        for entry in self.list_ensembles:
            self.system.ensembles.append(entry)

    def add_flr_fps_global_parameters(self,
                                      sheet_name='FLR_FPS_global_parameters',
                                      skiprows=3,
                                      header=0):
        """Read the flrCIF global parameters for fps."""
        if self.verbose:
            print(' ... Processing tab \'FLR_FPS_global_parameters\' ...')
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['FLR_FPS_global_parameters_Ordinal_id'])

        for i in range(nr_entries):
            cur_fps_glob_param_ordinal_id = \
                xls_data['FLR_FPS_global_parameters_Ordinal_id'][i]
            cur_fps_glob_param_ihm_modeling_protocol_id = \
                xls_data['FLR_FPS_global_parameters_IHM_modeling_protocol_id'][i]
            cur_fps_glob_param_ihm_modeling_protocol_step_id = \
                xls_data['FLR_FPS_global_parameters_IHM_modeling_protocol_step_id'][i]
            cur_fps_glob_param_fret_forster_radius_id = \
                xls_data['FLR_FPS_global_parameters_FRET_forster_radius'][i]
            cur_fps_glob_param_conversion_function_polynom_order = \
                xls_data['FLR_FPS_global_parameters_Conversion_function_polynom_order'][i]
            cur_fps_glob_param_repetition = \
                xls_data['FLR_FPS_global_parameters_Repetition'][i]
            cur_fps_glob_param_av_grid_rel = \
                xls_data['FLR_FPS_global_parameters_AV_grid_rel'][i]
            cur_fps_glob_param_av_min_grid_a = \
                xls_data['FLR_FPS_global_parameters_AV_min_grid_A'][i]
            cur_fps_glob_param_av_allowed_sphere = \
                xls_data['FLR_FPS_global_parameters_AV_allowed_sphere'][i]
            cur_fps_glob_param_av_search_nodes = \
                xls_data['FLR_FPS_global_parameters_AV_search_nodes'][i]
            cur_fps_glob_param_av_e_samples_k = \
                xls_data['FLR_FPS_global_parameters_AV_E_samples_k'][i]
            cur_fps_glob_param_sim_viscosity_adjustment = \
                None if ('FLR_FPS_global_parameters_Sim_viscosity_adjustment' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_viscosity_adjustment'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_viscosity_adjustment'][i]
            cur_fps_glob_param_sim_dt_adjustment = \
                None if ('FLR_FPS_global_parameters_Sim_dt_adjustment' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_dt_adjustment'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_dt_adjustment'][i]
            cur_fps_glob_param_sim_max_iter_k = \
                None if ('FLR_FPS_global_parameters_Sim_max_iter_k' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_max_iter_k'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_max_iter_k'][i]
            cur_fps_glob_param_sim_max_force = \
                None if ('FLR_FPS_global_parameters_Sim_max_force' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_max_force'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_max_force'][i]
            cur_fps_glob_param_sim_clash_tolerance_a = \
                None if ('FLR_FPS_global_parameters_Sim_clash_tolerance_A' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_clash_tolerance_A'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_clash_tolerance_A'][i]
            cur_fps_glob_param_sim_reciprocal_kt = \
                None if ('FLR_FPS_global_parameters_Sim_reciprocal_kT' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_reciprocal_kT'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_reciprocal_kT'][i]
            cur_fps_glob_param_sim_clash_potential = \
                None if ('FLR_FPS_global_parameters_Sim_clash_potential' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_Sim_clash_potential'][i])) \
                    else xls_data['FLR_FPS_global_parameters_Sim_clash_potential'][i]
            cur_fps_glob_param_convergence_e = \
                None if ('FLR_FPS_global_parameters_convergence_E' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_convergence_E'][i])) \
                    else xls_data['FLR_FPS_global_parameters_convergence_E'][i]
            cur_fps_glob_param_convergence_k = \
                None if ('FLR_FPS_global_parameters_convergence_K' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_convergence_K'][i])) \
                    else xls_data['FLR_FPS_global_parameters_convergence_K'][i]
            cur_fps_glob_param_convergence_f = \
                None if ('FLR_FPS_global_parameters_convergence_F' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_convergence_F'][i])) \
                    else xls_data['FLR_FPS_global_parameters_convergence_F'][i]
            cur_fps_glob_param_convergence_t = \
                None if ('FLR_FPS_global_parameters_convergence_T' not in
                         xls_data.keys() or pandas.isnull(
                            xls_data['FLR_FPS_global_parameters_convergence_T'][i])) \
                    else xls_data['FLR_FPS_global_parameters_convergence_T'][i]

            cur_fps_glob_param = ihm.flr.FPSGlobalParameters(
                forster_radius=
                cur_fps_glob_param_fret_forster_radius_id,
                conversion_function_polynom_order=
                cur_fps_glob_param_conversion_function_polynom_order,
                repetition=cur_fps_glob_param_repetition,
                av_grid_rel=cur_fps_glob_param_av_grid_rel,
                av_min_grid_a=cur_fps_glob_param_av_min_grid_a,
                av_allowed_sphere=cur_fps_glob_param_av_allowed_sphere,
                av_search_nodes=cur_fps_glob_param_av_search_nodes,
                av_e_samples_k=cur_fps_glob_param_av_e_samples_k,
                sim_viscosity_adjustment=
                cur_fps_glob_param_sim_viscosity_adjustment,
                sim_dt_adjustment=cur_fps_glob_param_sim_dt_adjustment,
                sim_max_iter_k=cur_fps_glob_param_sim_max_iter_k,
                sim_max_force=cur_fps_glob_param_sim_max_force,
                sim_clash_tolerance_a=cur_fps_glob_param_sim_clash_tolerance_a,
                sim_reciprocal_kt=cur_fps_glob_param_sim_reciprocal_kt,
                sim_clash_potential=cur_fps_glob_param_sim_clash_potential,
                convergence_e=cur_fps_glob_param_convergence_e,
                convergence_f=cur_fps_glob_param_convergence_f,
                convergence_k=cur_fps_glob_param_convergence_k,
                convergence_t=cur_fps_glob_param_convergence_t)

            if cur_fps_glob_param not in self.list_flr_fps_global_parameters:
                self.list_flr_fps_global_parameters.append(cur_fps_glob_param)
            # store the fps global parameters by the protocol id and the step id
            if cur_fps_glob_param_ihm_modeling_protocol_id not in \
                    self.dict_flr_fps_global_parameters_by_protocol_id.keys():
                self.dict_flr_fps_global_parameters_by_protocol_id[
                    cur_fps_glob_param_ihm_modeling_protocol_id] = {}
            if cur_fps_glob_param_ihm_modeling_protocol_step_id not in \
                    self.dict_flr_fps_global_parameters_by_protocol_id[
                        cur_fps_glob_param_ihm_modeling_protocol_id].keys():
                self.dict_flr_fps_global_parameters_by_protocol_id[
                    cur_fps_glob_param_ihm_modeling_protocol_id][
                    cur_fps_glob_param_ihm_modeling_protocol_step_id] = None
            # add the global parameters to the list by the key of the protocol
            # and the protocol step
            self.dict_flr_fps_global_parameters_by_protocol_id[
                cur_fps_glob_param_ihm_modeling_protocol_id][
                cur_fps_glob_param_ihm_modeling_protocol_step_id] = \
                cur_fps_glob_param
    def add_flr_fps_mpp(self,
                        sheet_name='FLR_FPS_MPP_group',
                        skiprows=3,
                        header=0):
        """Read the mean position for FPS from the excel sheet.
         In FPS it is possible to define mean positions for the AV center
         instead of creating the AVs.
         This is not recommended, but possible.
         The groups are stored as a list of atoms in a dictionary using the
         group id as key."""
        # First check whether the sheet with name 'FLR_FPS_MPP_group' exists
        if sheet_name in self.xls_file.sheet_names:
            if self.verbose:
                print(' ... Processing tab \'FPS_FPS_MPP_group\' ...')
            xls_data = pandas.read_excel(self.xls_file,
                                         sheet_name=sheet_name,
                                         skiprows=skiprows,
                                         header=header)

            nr_entries = len(xls_data['FPS_MPP_group_Ordinal_id'])

            for i in range(nr_entries):
                cur_fps_mpp_group_ordinal_id = \
                    xls_data['FPS_MPP_group_Ordinal_id'][i]
                cur_fps_mpp_group_id = \
                    xls_data['FPS_MPP_group_ID'][i]
                cur_fps_mpp_group_atom_position_entity_id = \
                    xls_data['FPS_MPP_atom_position_Entity'][i]
                cur_fps_mpp_group_atom_position_seq_id = \
                    xls_data['FPS_MPP_atom_position_Seq_id'][i]
                cur_fps_mpp_group_atom_position_chem_comp_id = \
                    xls_data['FPS_MPP_atom_position_Chem_comp_id'][i]
                cur_fps_mpp_atom_position_atom_id = \
                    xls_data['FPS_MPP_atom_position_Atom_id'][i]
                cur_fps_mpp_atom_position_asym_id = \
                    xls_data['FPS_MPP_atom_position_Asym_id'][i]
                cur_fps_mpp_atom_position_xcoord = \
                    xls_data['FPS_MPP_atom_position_Xcoord'][i]
                cur_fps_mpp_atom_position_ycoord = \
                    xls_data['FPS_MPP_atom_position_Ycoord'][i]
                cur_fps_mpp_atom_position_zcoord = \
                    xls_data['FPS_MPP_atom_position_Zcoord'][i]

                # Create the Residue or atom
                cur_resatom_new = None
                cur_entity = self.list_ihm_entities[
                    self.list_ihm_entity_ids.index(
                        cur_fps_mpp_group_atom_position_entity_id)]
                ## First create the residue
                cur_asym_unit = self.list_asym_units[
                    self.list_asym_units_ids.index(
                        cur_fps_mpp_atom_position_asym_id)]
                cur_residue = cur_asym_unit.residue(
                    seq_id=cur_fps_mpp_group_atom_position_seq_id)

                cur_atom = cur_residue.atom(
                    atom_id=cur_fps_mpp_atom_position_atom_id)
                cur_resatom_new = cur_atom
                # add it to the list of resatoms if it is not there yet
                if get_resatom_from_list(cur_resatom_new, self.list_resatoms) is None:
                    self.list_resatoms.append(cur_resatom_new)
                # and get the respective entry (to avoid duplicate entries)
                cur_resatom = get_resatom_from_list(cur_resatom_new, self.list_resatoms)

                cur_fps_mpp_atom_position = ihm.flr.FPSMPPAtomPosition(
                    atom=cur_resatom,
                    x=cur_fps_mpp_atom_position_xcoord,
                    y=cur_fps_mpp_atom_position_ycoord,
                    z=cur_fps_mpp_atom_position_zcoord)
                # if the current group is not in the list of mpp groups yet, add it
                if not cur_fps_mpp_group_id in self.dict_flr_fps_mpp_groups.keys():
                    self.dict_flr_fps_mpp_groups[cur_fps_mpp_group_id] = \
                        ihm.flr.FPSMPPAtomPositionGroup()
                self.dict_flr_fps_mpp_groups[
                    cur_fps_mpp_group_id].add_atom_position(
                    cur_fps_mpp_atom_position)

    def add_flr_reference_measurements(self,
                                       sheet_name='Reference_measurements',
                                       skiprows=3,
                                       header=0):
        """Read information on the reference measurements from the excel sheet"""
        if self.verbose:
            print(' ... Processing tab \'Reference_measurements\' ...')
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['RefMeas_ID'])

        for i in range(nr_entries):
            self.list_of_object_indices_refmeas.append({})

        ## only do the rest if there are reference measurements
        if nr_entries == 0:
            return
        # List that contains, whether the reference measurement contains
        # information about donor or acceptor
        tmp_list_ref_measurement_contains_donor_info = []
        tmp_list_ref_measurement_contains_acceptor_info = []
        for i in range(nr_entries):
            if not pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_donor_seq_id'][i]):
                tmp_list_ref_measurement_contains_donor_info.append(True)
            else:
                tmp_list_ref_measurement_contains_donor_info.append(False)
            if not pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_acceptor_seq_id'][i]):
                tmp_list_ref_measurement_contains_acceptor_info.append(True)
            else:
                tmp_list_ref_measurement_contains_acceptor_info.append(False)

        # Ensure that only either donor or acceptor are used
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i] and \
                    tmp_list_ref_measurement_contains_acceptor_info[i]:
                raise ValueError(
                    "ERROR: Reference measurements should only be used for "
                    "either donor or acceptor. If a sample was used for "
                    "both, please add a separate entry in the"
                    " \'Reference_measurements\' tab.")

        # Create the objects
        # => Sample conditions for reference measurements
        for i in range(nr_entries):
            # get the respective column for the sample condition
            cur_sample_condition_details = \
                None if ('RefMeas_Sample_condition' not in xls_data.keys() or
                         pandas.isnull(
                             xls_data['RefMeas_Sample_condition'][i])) \
                    else xls_data['RefMeas_Sample_condition'][i]
            # create the sample condition object
            cur_sample_condition = ihm.flr.SampleCondition(
                details=cur_sample_condition_details)
            # check whether it is already in the list
            cur_sample_condition_index = -1
            if cur_sample_condition not in self.list_sample_conditions:
                self.list_sample_conditions.append(cur_sample_condition)
            # and store the index of the respective object
            cur_sample_condition_index = \
                self.list_sample_conditions.index(cur_sample_condition)
            self.list_of_object_indices_refmeas[i]['Sample_condition'] = \
                cur_sample_condition_index

#            # => Instrument for reference measurements
#            for i in range(nr_entries):
#                # get the respective column for the sample condition
#                cur_instrument_details = \
#                    None if ('RefMeas_Instrument' not in  xls_data.keys() or
#                             pandas.isnull(xls_data['RefMeas_Instrument'][i])) \
#                        else xls_data['RefMeas_Instrument'][i]
#                ## create the object
#                cur_instrument = ihm.flr.Instrument(details=cur_instrument_details)
#                ## check whether it is already in the list
#                cur_instrument_index = -1
#                if cur_instrument not in self.list_instruments:
#                    self.list_instruments.append(cur_instrument)
#                ## and store the index of the respective object
#                cur_instrument_index = self.list_instruments.index(cur_instrument)
#                self.list_of_object_indices_refmeas[i]['Instrument'] = cur_instrument_index
#
#            # => Instrument settings
#            for i in range(nr_entries):
#                cur_inst_setting_details = \
#                    None if ('RefMeas_Instrument_setting' not in
#                             xls_data.keys() or pandas.isnull(
#                                xls_data['RefMeas_Instrument_setting'][i])) \
#                        else xls_data['RefMeas_Instrument_setting'][i]
#                cur_inst_setting = ihm.flr.InstSetting(details=cur_inst_setting_details)
#                cur_inst_setting_index = -1
#                if cur_inst_setting not in self.list_inst_settings:
#                    self.list_inst_settings.append(cur_inst_setting)
#                cur_inst_setting_index = self.list_inst_settings.index(cur_inst_setting)
#                self.list_of_object_indices_refmeas[i]['Instrument_setting'] = \
#                    cur_inst_setting_index
#
#            # => Experimental condition
#            for i in range(nr_entries):
#                cur_exp_condition_details = \
#                    None if ('RefMeas_Experimental_condition' not in
#                             xls_data.keys() or pandas.isnull(
#                                xls_data['RefMeas_Experimental_condition'][i])) \
#                        else xls_data['RefMeas_Experimental_condition'][i]
#                cur_exp_condition = ihm.flr.ExpCondition(details=cur_exp_condition_details)
#                cur_exp_condition_index = -1
#                if cur_exp_condition not in self.list_exp_conditions:
#                    self.list_exp_conditions.append(cur_exp_condition)
#                cur_exp_condition_index = self.list_exp_conditions.index(cur_exp_condition)
#                self.list_of_object_indices_refmeas[i][
#                    'Experimental_condition'] = cur_exp_condition_index
#

        # => Samples for reference measurements
        for i in range(nr_entries):
            cur_sample_id = xls_data['RefMeas_Sample_ID']
            cur_sample_num_of_probes = \
                int(xls_data['RefMeas_Sample_num_of_probes'][i])
            cur_sample_solvent_phase = \
                xls_data['RefMeas_Sample_solvent_phase'][i]
            cur_sample_description = \
                None if ('RefMeas_Sample_description' not in xls_data.keys()
                         or pandas.isnull(xls_data['RefMeas_Sample_description'][i])) \
                    else xls_data['RefMeas_Sample_description'][i]
            cur_sample_details = \
                None if ('RefMeas_Sample_details' not in xls_data.keys()
                         or pandas.isnull(xls_data['RefMeas_Sample_details'][i])) \
                    else xls_data['RefMeas_Sample_details'][i]
            cur_entity_assembly = xls_data['RefMeas_Entity_assembly'][i]
            # create the object
            # The entity assembly is taken from the list of entity assemblies
            # by getting the entry that corresponds to the given
            # entity_assembly_id (cur_entity_assembly)
            cur_sample = ihm.flr.Sample(
                entity_assembly=self.list_entity_assemblies[
                    self.list_entity_assembly_ids.index(cur_entity_assembly)],
                num_of_probes=cur_sample_num_of_probes,
                condition=self.list_sample_conditions[
                    self.list_of_object_indices_refmeas[i]['Sample_condition']],
                description=cur_sample_description,
                details=cur_sample_details,
                solvent_phase=cur_sample_solvent_phase)
            if not occurs_in_list(cur_sample, self.list_samples):
                self.list_samples.append(cur_sample)
                self.list_sample_ids.append(cur_sample_id)

            cur_sample_index = self.list_samples.index(cur_sample)
            self.list_of_object_indices_refmeas[i]['Sample'] = cur_sample_index

        # => Experiment for reference measurement
#        ## in case of the experiment, there is no need to read something from the file. Everything is defined before
#        Experiment_1 = ihm.flr.Experiment()
#        cur_experiment_index = -1
#        for i in range(nr_of_entries_refmeas):
#            ## make sure that every instrument-inst_setting-exp_condition-sample-combination is there only once
#            cur_instrument = list_instruments[list_of_object_indices_refmeas[i]['Instrument']]
#            cur_inst_setting = list_inst_settings[list_of_object_indices_refmeas[i]['Instrument_setting']]
#            cur_exp_condition = list_exp_conditions[list_of_object_indices_refmeas[i]['Experimental_condition']]
#            cur_sample = list_samples[list_of_object_indices_refmeas[i]['Sample']]
#            cur_experiment_details = None if (
#                   'RefMeas_Experiment_Details' not in xls_refmeas_data.keys() or pandas.isnull(
#                xls_refmeas_data['RefMeas_Experiment_Details'][i]))
#                       else xls_refmeas_data['RefMeas_Experiment_Details'][i]
#
#            if not Experiment_1.contains(cur_instrument, cur_inst_setting, cur_exp_condition, cur_sample):
#                Experiment_1.add_entry(instrument=cur_instrument, inst_setting=cur_inst_setting,
#                                       exp_condition=cur_exp_condition, sample=cur_sample,
#                                       details=cur_experiment_details)
#        if Experiment_1 not in list_experiments:
#            list_experiments.append(Experiment_1)
#        for i in range(nr_of_entries_refmeas):
#            cur_experiment_index = list_experiments.index(Experiment_1)
#            list_of_object_indices_refmeas[i]['Experiment'] = cur_experiment_index

        # => Probe for reference measurement

        # => Donor
        # Probe_list
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i] == True:
                cur_probe_donor_name = \
                    xls_data['RefMeas_Probe_donor_name'][i]
                cur_probe_donor_origin = \
                    xls_data['RefMeas_Probe_donor_origin'][i]
                cur_probe_donor_link_type = \
                    xls_data['RefMeas_Probe_donor_link_type'][i]
                cur_reactive_probe_donor_flag = False
                cur_reactive_probe_donor_name = None
                ## TODO! MODIFY TO HANDLE EMPTY CELLS HERE!
#                if 'RefMeas_Probe_donor_name' in xls_refmeas_data.keys():
                if 'RefMeas_Reactive_probe_donor_name' in xls_data.keys() \
                        and not pandas.isnull(xls_data['RefMeas_Probe_donor_name'][i]):
                    cur_reactive_probe_donor_flag = True
                    cur_reactive_probe_donor_name = \
                        xls_data['RefMeas_Reactive_probe_donor_name'][i]
                cur_probe_list_donor = \
                    ihm.flr.ProbeList(
                        chromophore_name=cur_probe_donor_name,
                        reactive_probe_flag=cur_reactive_probe_donor_flag,
                        reactive_probe_name=cur_reactive_probe_donor_name,
                        probe_origin=cur_probe_donor_origin,
                        probe_link_type=cur_probe_donor_link_type)
                # Probe_descriptor
                cur_probe_descriptor_donor = None
                cur_reactive_probe_donor_chemical_descriptor = None
                # index of the reactive_probe_chemical_descriptor in the
                # list_chemical_descriptors
                cur_chromophore_donor_chemical_descriptor = None
                # index of the chromophore in the list_chemical_descriptors
                # Reactive probe donor chemical descriptor
                # if there is a name given and at least a smiles,
                # canonical smiles, inchi, or inchi-key
                if ('RefMeas_Reactive_probe_donor_name' in xls_data.keys()
                        and not pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_name'][i])
                        and (('RefMeas_Reactive_probe_donor_smiles' in xls_data.keys()
                              and not pandas.isnull(
                                    xls_data['RefMeas_Reactive_probe_donor_smiles'][i]))
                             or ('RefMeas_Reactive_probe_donor_smiles_canonical' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Reactive_probe_donor_smiles_canonical'][i]))
                             or ('RefMeas_Reactive_probe_donor_inchi' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Reactive_probe_donor_inchi'][i]))
                             or ('RefMeas_Reactive_probe_donor_inchi_key' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Reactive_probe_donor_inchi_key'][i])))):
                    cur_reactive_probe_donor_chemical_descriptor = ihm.ChemDescriptor(
                        auth_name=xls_data['RefMeas_Reactive_probe_donor_name'][i],
                        chem_comp_id=None if (
                                'RefMeas_Reactive_probe_donor_chem_comp_id' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_chem_comp_id'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_chem_comp_id'][i],
                        chemical_name=None if (
                                'RefMeas_Reactive_probe_donor_chemical_name' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_chemical_name'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_chemical_name'][i],
                        common_name=None if (
                                'RefMeas_Reactive_probe_donor_common_name' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_common_name'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_common_name'][i],
                        smiles=None if (
                                'RefMeas_Reactive_probe_donor_smiles' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_smiles'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_smiles'][i],
                        smiles_canonical=None if (
                                'RefMeas_Reactive_probe_donor_smiles_canonical' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_smiles_canonical'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_smiles_canonical'][i],
                        inchi=None if (
                                'RefMeas_Reactive_probe_donor_inchi' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_donor_inchi'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_inchi'][i],
                        inchi_key=None if (
                                'RefMeas_Reactive_probe_donor_inchi_key' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_donor_inchi_key'][i]))
                        else xls_data['RefMeas_Reactive_probe_donor_inchi_key'][i])

                # Same for the chromophore - name and some structural
                # description (smiles or canonical smiles or inchi or inchi key)
                if ('RefMeas_Chromophore_donor_name' in xls_data.keys()
                        and not pandas.isnull(xls_data['RefMeas_Chromophore_donor_name'][i])
                        and (('RefMeas_Chromophore_donor_smiles' in xls_data.keys()
                              and not pandas.isnull(
                                    xls_data['RefMeas_Chromophore_donor_smiles'][i]))
                             or ('RefMeas_Chromophore_donor_smiles_canonical' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Chromophore_donor_smiles_canonical'][i]))
                             or ('RefMeas_Chromophore_donor_inchi' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Chromophore_donor_inchi'][i]))
                             or ('RefMeas_Chromophore_donor_inchi_key' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Chromophore_donor_inchi_key'][i])))):
                    cur_chromophore_donor_chemical_descriptor = \
                        ihm.ChemDescriptor(
                            auth_name=xls_data['RefMeas_Chromophore_donor_name'][i],
                            chem_comp_id=None if (
                                    'RefMeas_Chromophore_donor_chem_comp_id' not in
                                    xls_data.keys() or pandas.isnull(
                                xls_data['RefMeas_Chromophore_donor_chem_comp_id'][i]))
                            else xls_data['RefMeas_Chromophore_donor_chem_comp_id'][i],
                            chemical_name=None if (
                                    'RefMeas_Chromophore_donor_chemical_name' not in xls_data.keys()
                                    or pandas.isnull(
                                xls_data['RefMeas_Chromophore_donor_chemical_name'][i]))
                            else xls_data['RefMeas_Chromophore_donor_chemical_name'][i],
                            common_name=None if (
                                    'RefMeas_Chromophore_donor_common_name' not in xls_data.keys()
                                    or pandas.isnull(
                                xls_data['RefMeas_Chromophore_donor_common_name'][i]))
                            else xls_data['RefMeas_Chromophore_donor_common_name'][i],
                            smiles=None if (
                                    'RefMeas_Chromophore_donor_smiles' not in xls_data.keys()
                                    or pandas.isnull(
                                xls_data['RefMeas_Chromophore_donor_smiles'][i]))
                            else xls_data['RefMeas_Chromophore_donor_smiles'][i],
                            smiles_canonical=None if (
                                    'RefMeas_Chromophore_donor_smiles_canonical' not in xls_data.keys()
                                    or pandas.isnull(
                                xls_data['RefMeas_Chromophore_donor_smiles_canonical'][i]))
                            else xls_data['RefMeas_Chromophore_donor_smiles_canonical'][i],
                            inchi=None if (
                                    'RefMeas_Chromophore_donor_inchi' not in xls_data.keys()
                                    or pandas.isnull(
                                xls_data['RefMeas_Chromophore_donor_inchi'][i]))
                            else xls_data['RefMeas_Chromophore_donor_inchi'][i],
                            inchi_key=None if (
                                    'RefMeas_Chromophore_donor_inchi_key' not in xls_data.keys()
                                    or pandas.isnull(xls_data['RefMeas_Chromophore_donor_inchi_key'][i]))
                            else xls_data['RefMeas_Chromophore_donor_inchi_key'][i])

                # In case the current chemical descriptor is already present
                # in the chemical descriptors list, we use that one instead of
                # creating a new one.
                if cur_reactive_probe_donor_chemical_descriptor is not None:
                    if not occurs_in_list(
                            cur_reactive_probe_donor_chemical_descriptor,
                            self.list_chemical_descriptors):
                        self.list_chemical_descriptors.append(
                            cur_reactive_probe_donor_chemical_descriptor)
                    else:
                        cur_reactive_probe_donor_chemical_descriptor = \
                            [x for x in self.list_chemical_descriptors if
                             x.__dict__ == cur_reactive_probe_donor_chemical_descriptor.__dict__][0]
                if cur_chromophore_donor_chemical_descriptor is not None:
                    if not occurs_in_list(
                            cur_chromophore_donor_chemical_descriptor,
                            self.list_chemical_descriptors):
                        self.list_chemical_descriptors.append(
                            cur_chromophore_donor_chemical_descriptor)
                    else:
                        cur_chromophore_donor_chemical_descriptor = \
                            [x for x in self.list_chemical_descriptors if
                             x.__dict__ == cur_chromophore_donor_chemical_descriptor.__dict__][0]

                cur_probe_descriptor_donor = ihm.flr.ProbeDescriptor(
                    reactive_probe_chem_descriptor=
                    cur_reactive_probe_donor_chemical_descriptor,
                    chromophore_chem_descriptor=
                    cur_chromophore_donor_chemical_descriptor,
                    chromophore_center_atom=
                    None if (
                            'RefMeas_Chromophore_donor_center_atom' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Chromophore_donor_center_atom'][i]))
                    else xls_data['RefMeas_Chromophore_donor_center_atom'][i])

                if not occurs_in_list(
                        cur_probe_descriptor_donor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(
                        cur_probe_descriptor_donor)
                else:
                    cur_probe_descriptor_donor = \
                        [x for x in self.list_chemical_descriptors if
                         x.__dict__ == cur_probe_descriptor_donor.__dict__][0]

                cur_probe_donor = ihm.flr.Probe(
                    probe_list_entry=cur_probe_list_donor,
                    probe_descriptor=cur_probe_descriptor_donor)

                cur_probe_donor_index = -1
                if cur_probe_donor not in self.list_probe_donors:
                    self.list_probe_donors.append(cur_probe_donor)
                cur_probe_donor_index = self.list_probe_donors.index(cur_probe_donor)
                self.list_of_object_indices_refmeas[i]['Probe_donor'] = \
                    cur_probe_donor_index

        # => Acceptor
        # Probe_list
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i] == True:
                cur_probe_acceptor_name = \
                    xls_data['RefMeas_Probe_acceptor_name'][i]
                cur_probe_acceptor_origin = \
                    xls_data['RefMeas_Probe_acceptor_origin'][i]
                cur_probe_acceptor_link_type = \
                    xls_data['RefMeas_Probe_acceptor_link_type'][i]
                cur_reactive_probe_acceptor_flag = False
                cur_reactive_probe_acceptor_name = None
                if 'RefMeas_Reactive_probe_acceptor_name' in xls_data.keys() \
                        and not pandas.isnull(
                    xls_data['RefMeas_Reactive_probe_acceptor_name'][i]):
                    cur_reactive_probe_acceptor_flag = True
                    cur_reactive_probe_acceptor_name = \
                        xls_data['RefMeas_Reactive_probe_acceptor_name'][i]
                cur_probe_list_acceptor = \
                    ihm.flr.ProbeList(
                        chromophore_name=cur_probe_acceptor_name,
                        reactive_probe_flag=cur_reactive_probe_acceptor_flag,
                        reactive_probe_name=cur_reactive_probe_acceptor_name,
                        probe_origin=cur_probe_acceptor_origin,
                        probe_link_type=cur_probe_acceptor_link_type)
                # Probe_descriptor
                cur_probe_descriptor_acceptor = None
                cur_reactive_probe_acceptor_chemical_descriptor = None
                # index of the reactive_probe_chemical_descriptor in the
                # list_chemical_descriptors
                cur_chromophore_acceptor_chemical_descriptor = None
                # Reactive probe acceptor chemical descriptor
                # if there is a name given and at least a smiles,
                # canonical smiles, inchi, or inchi-key
                if ('RefMeas_Reactive_probe_acceptor_name' in xls_data.keys()
                        and not pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_acceptor_name'][i])
                        and (('RefMeas_Reactive_probe_acceptor_smiles' in xls_data.keys()
                              and not pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_smiles'][i]))
                             or ('RefMeas_Reactive_probe_acceptor_smiles_canonical' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Reactive_probe_acceptor_smiles_canonical'][i]))
                             or ('RefMeas_Reactive_probe_acceptor_inchi' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Reactive_probe_acceptor_inchi'][i]))
                             or ('RefMeas_Reactive_probe_acceptor_inchi_key' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Reactive_probe_acceptor_inchi_key'][i])))):
                    cur_reactive_probe_acceptor_chemical_descriptor = ihm.ChemDescriptor(
                        auth_name=xls_data['RefMeas_Reactive_probe_acceptor_name'][i],
                        chem_comp_id=None if (
                                'RefMeas_Reactive_probe_acceptor_chem_comp_id' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_chem_comp_id'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_chem_comp_id'][i],
                        chemical_name=None if (
                                'RefMeas_Reactive_probe_acceptor_chemical_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_chemical_name'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_chemical_name'][i],
                        common_name=None if (
                                'RefMeas_Reactive_probe_acceptor_common_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_common_name'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_common_name'][i],
                        smiles=None if (
                                'RefMeas_Reactive_probe_acceptor_smiles' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_smiles'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_smiles'][i],
                        smiles_canonical=None if (
                                'RefMeas_Reactive_probe_acceptor_smiles_canonical' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Reactive_probe_acceptor_smiles_canonical'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_smiles_canonical'][i],
                        inchi=None if (
                                'RefMeas_Reactive_probe_acceptor_inchi' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_inchi'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_inchi'][i],
                        inchi_key=None if (
                                'RefMeas_Reactive_probe_acceptor_inchi_key' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Reactive_probe_acceptor_inchi_key'][i]))
                        else xls_data['RefMeas_Reactive_probe_acceptor_inchi_key'][i])

                if ('RefMeas_Chromophore_acceptor_name' in xls_data.keys()
                        and not pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_name'][i])
                        and (('RefMeas_Chromophore_acceptor_smiles' in xls_data.keys()
                              and not pandas.isnull(
                                    xls_data['RefMeas_Chromophore_acceptor_smiles'][i]))
                             or ('RefMeas_Chromophore_acceptor_smiles_canonical' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Chromophore_acceptor_smiles_canonical'][i]))
                             or ('RefMeas_Chromophore_acceptor_inchi' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Chromophore_acceptor_inchi'][i]))
                             or ('RefMeas_Chromophore_acceptor_inchi_key' in xls_data.keys()
                                 and not pandas.isnull(
                                            xls_data['RefMeas_Chromophore_acceptor_inchi_key'][i])))):
                    cur_chromophore_acceptor_chemical_descriptor = ihm.ChemDescriptor(
                        auth_name=xls_data['RefMeas_Chromophore_acceptor_name'][i],
                        chem_comp_id=None if (
                                'RefMeas_Chromophore_acceptor_chem_comp_id' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_chem_comp_id'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_chem_comp_id'][i],
                        chemical_name=None if (
                                'RefMeas_Chromophore_acceptor_chemical_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_chemical_name'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_chemical_name'][i],
                        common_name=None if (
                                'RefMeas_Chromophore_acceptor_common_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_common_name'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_common_name'][i],
                        smiles=None if (
                                'RefMeas_Chromophore_acceptor_smiles' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_smiles'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_smiles'][i],
                        smiles_canonical=None if (
                                'RefMeas_Chromophore_acceptor_smiles_canonical' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_smiles_canonical'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_smiles_canonical'][i],
                        inchi=None if (
                                'RefMeas_Chromophore_acceptor_inchi' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_inchi'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_inchi'][i],
                        inchi_key=None if (
                                'RefMeas_Chromophore_acceptor_inchi_key' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_inchi_key'][i]))
                        else xls_data['RefMeas_Chromophore_acceptor_inchi_key'][i])

                if cur_reactive_probe_acceptor_chemical_descriptor is not None:
                    if not occurs_in_list(
                            cur_reactive_probe_acceptor_chemical_descriptor,
                            self.list_chemical_descriptors):
                        self.list_chemical_descriptors.append(
                            cur_reactive_probe_acceptor_chemical_descriptor)
                    else:
                        cur_reactive_probe_acceptor_chemical_descriptor = \
                            [x for x in self.list_chemical_descriptors if
                             x.__dict__ == cur_reactive_probe_acceptor_chemical_descriptor.__dict__][0]
                if cur_chromophore_acceptor_chemical_descriptor is not None:
                    if not occurs_in_list(
                            cur_chromophore_acceptor_chemical_descriptor,
                            self.list_chemical_descriptors):
                        self.list_chemical_descriptors.append(
                            cur_chromophore_acceptor_chemical_descriptor)
                    else:
                        cur_chromophore_acceptor_chemical_descriptor = \
                            [x for x in self.list_chemical_descriptors if
                             x.__dict__ == cur_chromophore_acceptor_chemical_descriptor.__dict__][0]

                cur_probe_descriptor_acceptor = ihm.flr.ProbeDescriptor(
                    reactive_probe_chem_descriptor=cur_reactive_probe_acceptor_chemical_descriptor,
                    chromophore_chem_descriptor=cur_chromophore_acceptor_chemical_descriptor,
                    chromophore_center_atom=
                    None if (
                            'RefMeas_Chromophore_acceptor_center_atom' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Chromophore_acceptor_center_atom'][i]))
                    else xls_data['RefMeas_Chromophore_acceptor_center_atom'][i])

                if not occurs_in_list(cur_probe_descriptor_acceptor, self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(cur_probe_descriptor_acceptor)
                else:
                    cur_probe_descriptor_acceptor = \
                        [x for x in self.list_chemical_descriptors if
                         x.__dict__ == cur_probe_descriptor_acceptor.__dict__][0]

                cur_probe_acceptor = ihm.flr.Probe(probe_list_entry=cur_probe_list_acceptor,
                                                   probe_descriptor=cur_probe_descriptor_acceptor)

                cur_probe_acceptor_index = -1
                if cur_probe_acceptor not in self.list_probe_acceptors:
                    self.list_probe_acceptors.append(cur_probe_acceptor)
                cur_probe_acceptor_index = self.list_probe_acceptors.index(cur_probe_acceptor)
                self.list_of_object_indices_refmeas[i]['Probe_acceptor'] = cur_probe_acceptor_index

        # => Modified and mutated residues
        # donor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i]:
                if ('RefMeas_Modified_residue_donor_name' in xls_data.keys()
                    and not pandas.isnull(xls_data['RefMeas_Modified_residue_donor_name'][i])) and \
                        (('RefMeas_Modified_residue_donor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['RefMeas_Modified_residue_donor_smiles'][i]))
                         or ('RefMeas_Modified_residue_donor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(
                                            xls_data['RefMeas_Modified_residue_donor_smiles_canonical'][i]))
                         or ('RefMeas_Modified_residue_donor_inchi' in xls_data.keys()
                             and not pandas.isnull(xls_data['RefMeas_Modified_residue_donor_inchi'][i]))
                         or ('RefMeas_Modified_residue_donor_inchi_key' in xls_data.keys()
                             and not pandas.isnull(xls_data['RefMeas_Modified_residue_donor_inchi_key'][i]))):
                    cur_modified_residue_donor_chemical_descriptor = ihm.ChemDescriptor(
                        auth_name=xls_data['RefMeas_Modified_residue_donor_name'][i],
                        chem_comp_id=None if (
                                'RefMeas_Modified_residue_donor_chem_comp_id' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_donor_chem_comp_id'][i]))
                        else xls_refmeas_data['RefMeas_Modified_residue_donor_chem_comp_id'][i],
                        chemical_name=None if (
                                'RefMeas_Modified_residue_donor_chemical_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_donor_chemical_name'][i]))
                        else xls_data['RefMeas_Modified_residue_donor_chemical_name'][i],
                        common_name=None if (
                                'RefMeas_Modified_residue_donor_common_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_donor_common_name'][i]))
                        else xls_data['RefMeas_Modified_residue_donor_common_name'][i],
                        smiles=None if (
                                'RefMeas_Modified_residue_donor_smiles' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_donor_smiles'][i]))
                        else xls_data['RefMeas_Modified_residue_donor_smiles'][i],
                        smiles_canonical=None if (
                                'RefMeas_Modified_residue_donor_smiles_canonical' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Modified_residue_donor_smiles_canonical'][i]))
                        else xls_data['RefMeas_Modified_residue_donor_smiles_canonical'][i],
                        inchi=None if (
                                'RefMeas_Modified_residue_donor_inchi' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_donor_inchi'][i]))
                        else xls_data['RefMeas_Modified_residue_donor_inchi'][i],
                        inchi_key=None if (
                                'RefMeas_Modified_residue_donor_inchi_key' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_donor_inchi_key'][i]))
                        else xls_data['RefMeas_Modified_residue_donor_inchi_key'][i])

                    if not occurs_in_list(
                            cur_modified_residue_donor_chemical_descriptor,
                            self.list_chemical_descriptors):
                        self.list_chemical_descriptors.append(
                            cur_modified_residue_donor_chemical_descriptor)
                    else:
                        cur_modified_residue_donor_chemical_descriptor = [
                            x for x in self.list_chemical_descriptors
                            if x.__dict__ == cur_modified_residue_donor_chemical_descriptor.__dict__][0]

                    cur_modified_residue_donor_index = -1
                    if not occurs_in_list(
                            cur_modified_residue_donor_chemical_descriptor,
                            self.list_modified_residues):
                        self.list_modified_residues.append(
                            cur_modified_residue_donor_chemical_descriptor)
                    cur_modified_residue_donor_index = self.list_modified_residues.index(
                        cur_modified_residue_donor_chemical_descriptor)
                    self.list_of_object_indices_refmeas[i]['Modified_residue_donor'] = \
                        cur_modified_residue_donor_index

        # acceptor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                if ('RefMeas_Modified_residue_acceptor_name' in xls_data.keys()
                    and not pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_name'][i])) and \
                        (('RefMeas_Modified_residue_acceptor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_smiles'][i]))
                         or ('RefMeas_Modified_residue_acceptor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(
                                            xls_data['RefMeas_Modified_residue_acceptor_smiles_canonical'][i]))
                         or ('RefMeas_Modified_residue_acceptor_inchi' in xls_data.keys()
                             and not pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_inchi'][i]))
                         or ('RefMeas_Modified_residue_acceptor_inchi_key' in xls_data.keys()
                             and not pandas.isnull(
                                            xls_data['RefMeas_Modified_residue_acceptor_inchi_key'][i]))):
                    cur_modified_residue_acceptor_chemical_descriptor = ihm.ChemDescriptor(
                        auth_name=xls_data['RefMeas_Modified_residue_acceptor_name'][i],
                        chem_comp_id=None if (
                                'RefMeas_Modified_residue_acceptor_chem_comp_id' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_chem_comp_id'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_chem_comp_id'][i],
                        chemical_name=None if (
                                'RefMeas_Modified_residue_acceptor_chemical_name' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Modified_residue_acceptor_chemical_name'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_chemical_name'][i],
                        common_name=None if (
                                'RefMeas_Modified_residue_acceptor_common_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_common_name'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_common_name'][i],
                        smiles=None if (
                                'RefMeas_Modified_residue_acceptor_smiles' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_smiles'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_smiles'][i],
                        smiles_canonical=None if (
                                'RefMeas_Modified_residue_acceptor_smiles_canonical' not in xls_data.keys()
                                or pandas.isnull(
                            xls_data['RefMeas_Modified_residue_acceptor_smiles_canonical'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_smiles_canonical'][i],
                        inchi=None if (
                                'RefMeas_Modified_residue_acceptor_inchi' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_inchi'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_inchi'][i],
                        inchi_key=None if (
                                'RefMeas_Modified_residue_acceptor_inchi_key' not in xls_data.keys()
                                or pandas.isnull(xls_data['RefMeas_Modified_residue_acceptor_inchi_key'][i]))
                        else xls_data['RefMeas_Modified_residue_acceptor_inchi_key'][i])

                    if not occurs_in_list(
                            cur_modified_residue_acceptor_chemical_descriptor,
                            self.list_chemical_descriptors):
                        self.list_chemical_descriptors.append(
                            cur_modified_residue_acceptor_chemical_descriptor)
                    else:
                        cur_modified_residue_acceptor_chemical_descriptor = [
                            x for x in self.list_chemical_descriptors if
                            x.__dict__ == cur_modified_residue_acceptor_chemical_descriptor.__dict__][0]

                    cur_modified_residue_acceptor_index = -1
                    if not occurs_in_list(
                            cur_modified_residue_acceptor_chemical_descriptor,
                            self.list_modified_residues):
                        self.list_modified_residues.append(
                            cur_modified_residue_acceptor_chemical_descriptor)
                    cur_modified_residue_acceptor_index = self.list_modified_residues.index(
                        cur_modified_residue_acceptor_chemical_descriptor)
                    self.list_of_object_indices_refmeas[i]['Modified_residue_acceptor'] = \
                        cur_modified_residue_acceptor_index

        # ==> Mutated residues
        # donor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i]:
                if ('RefMeas_Mutated_residue_donor_ChemComp_ID' in xls_data.keys()
                        and not pandas.isnull(xls_data['RefMeas_Mutated_residue_donor_ChemComp_ID'][i])):
                    cur_mutated_residue_donor_chemcomp_id = \
                    xls_data['RefMeas_Mutated_residue_donor_ChemComp_ID'][i]
                    cur_mutated_residue_donor_code = None if (
                            'RefMeas_Mutated_residue_donor_Code' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Mutated_residue_donor_Code'][i])) \
                        else xls_data['RefMeas_Mutated_residue_donor_Code'][i]
                    cur_mutated_residue_donor_code_canonical = None if (
                            'RefMeas_Mutated_residue_donor_Code_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Mutated_residue_donor_Code_canonical'][i])) \
                        else xls_data['RefMeas_Mutated_residue_donor_Code_canonical'][i]
                    cur_mutated_residue_donor_type = None if (
                            'RefMeas_Mutated_residue_donor_Type' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Mutated_residue_donor_Type'][i])) \
                        else xls_data['RefMeas_Mutated_residue_donor_Type'][i]
                    cur_mutated_residue_donor_chem_comp = None
                    if cur_mutated_residue_donor_type == 'Peptide linking':
                        cur_mutated_residue_donor_chem_comp = ihm.PeptideChemComp(
                            id=cur_mutated_residue_donor_chemcomp_id,
                            code=cur_mutated_residue_donor_code,
                            code_canonical=cur_mutated_residue_donor_code_canonical)
                    elif cur_mutated_residue_donor_type == 'L-Peptide linking':
                        cur_mutated_residue_donor_chem_comp = ihm.LPeptideChemComp(
                            id=cur_mutated_residue_donor_chemcomp_id,
                            code=cur_mutated_residue_donor_code,
                            code_canonical=cur_mutated_residue_donor_code_canonical)
                    elif cur_mutated_residue_donor_type == 'D-Peptide linking':
                        cur_mutated_residue_donor_chem_comp = ihm.DPeptideChemComp(
                            id=cur_mutated_residue_donor_chemcomp_id,
                            code=cur_mutated_residue_donor_code,
                            code_canonical=cur_mutated_residue_donor_code_canonical)
                    elif cur_mutated_residue_donor_type == 'DNA linking':
                        cur_mutated_residue_donor_chem_comp = ihm.DNAChemComp(
                            id=cur_mutated_residue_donor_chemcomp_id,
                            code=cur_mutated_residue_donor_code,
                            code_canonical=cur_mutated_residue_donor_code_canonical)
                    elif cur_mutated_residue_donor_type == 'RNA linking':
                        cur_mutated_residue_donor_chem_comp = ihm.RNAChemComp(
                            id=cur_mutated_residue_donor_chemcomp_id,
                            code=cur_mutated_residue_donor_code,
                            code_canonical=cur_mutated_residue_donor_code_canonical)
                    else:
                        cur_mutated_residue_donor_chem_comp = ihm.ChemComp(
                            id=cur_mutated_residue_donor_chemcomp_id,
                            code=cur_mutated_residue_donor_code,
                            code_canonical=cur_mutated_residue_donor_code_canonical)

                    if not occurs_in_list(
                            cur_mutated_residue_donor_chem_comp,
                            self.list_chemical_components):
                        self.list_chemical_components.append(cur_mutated_residue_donor_chem_comp)
                    else:
                        cur_mutated_residue_donor_chem_comp = [
                            x for x in self.list_chemical_components if
                            x.__dict__ == cur_mutated_residue_donor_chem_comp.__dict__][0]

                    cur_mutated_residue_donor_index = -1
                    if not occurs_in_list(
                            cur_mutated_residue_donor_chem_comp,
                            self.list_mutated_residues):
                        self.list_mutated_residues.append(
                            cur_mutated_residue_donor_chem_comp)
                    cur_mutated_residue_donor_index = \
                        self.list_mutated_residues.index(cur_mutated_residue_donor_chem_comp)
                    self.list_of_object_indices_refmeas[i]['Mutated_residue_donor'] = \
                        cur_mutated_residue_donor_index

        # acceptor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                if ('RefMeas_Mutated_residue_acceptor_ChemComp_ID' in xls_data.keys()
                        and not pandas.isnull(xls_data['RefMeas_Mutated_residue_acceptor_ChemComp_ID'][i])):
                    cur_mutated_residue_acceptor_chemcomp_id = \
                    xls_data['RefMeas_Mutated_residue_acceptor_ChemComp_ID'][i]
                    cur_mutated_residue_acceptor_code = None if (
                            'RefMeas_Mutated_residue_acceptor_Code' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Mutated_residue_acceptor_Code'][i])) \
                        else xls_data['RefMeas_Mutated_residue_acceptor_Code'][i]
                    cur_mutated_residue_acceptor_code_canonical = None if (
                            'RefMeas_Mutated_residue_acceptor_Code_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Mutated_residue_acceptor_Code_canonical'][i])) \
                        else xls_data['RefMeas_Mutated_residue_acceptor_Code_canonical'][i]
                    cur_mutated_residue_acceptor_type = None if (
                            'RefMeas_Mutated_residue_acceptor_Type' not in xls_data.keys()
                            or pandas.isnull(xls_data['RefMeas_Mutated_residue_acceptor_Type'][i])) \
                        else xls_data['RefMeas_Mutated_residue_acceptor_Type'][i]
                    cur_mutated_residue_acceptor_chem_comp = None
                    if cur_mutated_residue_acceptor_type == 'Peptide linking':
                        cur_mutated_residue_acceptor_chem_comp = ihm.PeptideChemComp(
                            id=cur_mutated_residue_acceptor_chemcomp_id,
                            code=cur_mutated_residue_acceptor_code,
                            code_canonical=cur_mutated_residue_acceptor_code_canonical)
                    elif cur_mutated_residue_acceptor_type == 'L-Peptide linking':
                        cur_mutated_residue_acceptor_chem_comp = ihm.LPeptideChemComp(
                            id=cur_mutated_residue_acceptor_chemcomp_id,
                            code=cur_mutated_residue_acceptor_code,
                            code_canonical=cur_mutated_residue_acceptor_code_canonical)
                    elif cur_mutated_residue_acceptor_type == 'D-Peptide linking':
                        cur_mutated_residue_acceptor_chem_comp = ihm.DPeptideChemComp(
                            id=cur_mutated_residue_acceptor_chemcomp_id,
                            code=cur_mutated_residue_acceptor_code,
                            code_canonical=cur_mutated_residue_acceptor_code_canonical)
                    elif cur_mutated_residue_acceptor_type == 'DNA linking':
                        cur_mutated_residue_acceptor_chem_comp = ihm.DNAChemComp(
                            id=cur_mutated_residue_acceptor_chemcomp_id,
                            code=cur_mutated_residue_acceptor_code,
                            code_canonical=cur_mutated_residue_acceptor_code_canonical)
                    elif cur_mutated_residue_acceptor_type == 'RNA linking':
                        cur_mutated_residue_acceptor_chem_comp = ihm.RNAChemComp(
                            id=cur_mutated_residue_acceptor_chemcomp_id,
                            code=cur_mutated_residue_acceptor_code,
                            code_canonical=cur_mutated_residue_acceptor_code_canonical)
                    else:
                        cur_mutated_residue_acceptor_chem_comp = ihm.ChemComp(
                            id=cur_mutated_residue_acceptor_chemcomp_id,
                            code=cur_mutated_residue_acceptor_code,
                            code_canonical=cur_mutated_residue_acceptor_code_canonical)

                    if not occurs_in_list(
                            cur_mutated_residue_acceptor_chem_comp,
                            self.list_chemical_components):
                        self.list_chemical_components.append(
                            cur_mutated_residue_acceptor_chem_comp)
                    else:
                        cur_mutated_residue_acceptor_chem_comp = [
                            x for x in self.list_chemical_components if
                            x.__dict__ == cur_mutated_residue_acceptor_chem_comp.__dict__][0]

                    cur_mutated_residue_acceptor_index = -1
                    if not occurs_in_list(
                            cur_mutated_residue_acceptor_chem_comp,
                            self.list_mutated_residues):
                        self.list_mutated_residues.append(cur_mutated_residue_acceptor_chem_comp)
                    cur_mutated_residue_acceptor_index = self.list_mutated_residues.index(
                        cur_mutated_residue_acceptor_chem_comp)
                    self.list_of_object_indices_refmeas[i]['Mutated_residue_acceptor'] = \
                        cur_mutated_residue_acceptor_index

        # => Poly_probe_position
        # donor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i]:
                cur_poly_probe_position_auth_name = None if (
                        'RefMeas_Poly_probe_position_donor_name' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_donor_name'][i])) else \
                    xls_data['RefMeas_Poly_probe_position_donor_name'][i]
                cur_poly_probe_position_entity = \
                    xls_data['RefMeas_Poly_probe_position_donor_entity'][i]
                cur_poly_probe_position_asym_unit_id = None if (
                        'RefMeas_Poly_probe_position_donor_AsymUnit' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_donor_AsymUnit'][i])) else \
                    xls_data['RefMeas_Poly_probe_position_donor_AsymUnit'][i]
                cur_poly_probe_position_seq_id = \
                    int(xls_data['RefMeas_Poly_probe_position_donor_seq_id'][i])
                cur_poly_probe_position_comp_id = \
                    xls_data['RefMeas_Poly_probe_position_donor_comp_id'][i]
                cur_poly_probe_position_atom_id = None if (
                        'RefMeas_Poly_probe_position_donor_atom_id' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_donor_atom_id'][i])) else \
                    xls_data['RefMeas_Poly_probe_position_donor_atom_id'][i]
                cur_poly_probe_position_mutation_flag = False if (
                        'RefMeas_Poly_probe_position_donor_mutation_flag' not in xls_data.keys()
                        or xls_data['RefMeas_Poly_probe_position_donor_mutation_flag'][i] in
                        ['False', 'no', 'No', 'NO']) else True
                cur_poly_probe_position_modification_flag = False if (
                        'RefMeas_Poly_probe_position_donor_modification_flag' not in xls_data.keys()
                        or xls_data['RefMeas_Poly_probe_position_donor_modification_flag'][i] in
                        ['False', 'no', 'No', 'NO']) else True

                # Create the Residue or atom
                cur_resatom_new = None
                cur_entity = self.list_ihm_entities[
                    self.list_ihm_entity_ids.index(cur_poly_probe_position_entity)]
                cur_asym_unit = None if cur_poly_probe_position_asym_unit_id is None else \
                    self.list_asym_units[
                        self.list_asym_units_ids.index(cur_poly_probe_position_asym_unit_id)]
                # First create the residue
                cur_residue = cur_entity.residue(seq_id=cur_poly_probe_position_seq_id)
                if cur_asym_unit is not None:
                    cur_residue.asym = cur_asym_unit

                if cur_poly_probe_position_atom_id is not None:
                    cur_atom = cur_residue.atom(atom_id=cur_poly_probe_position_atom_id)
                    cur_resatom_new = cur_atom
                else:
                    cur_resatom_new = cur_residue
                # add it to the list of resatoms if it is not there yet
                if get_resatom_from_list(cur_resatom_new, self.list_resatoms) is None:
                    self.list_resatoms.append(cur_resatom_new)
                ## and get the respective entry (to avoid duplicate entries)
                cur_resatom = get_resatom_from_list(cur_resatom_new,
                                                    self.list_resatoms)

                cur_poly_probe_position = ihm.flr.PolyProbePosition(
                    resatom=cur_resatom,
                    mutation_flag=cur_poly_probe_position_mutation_flag,
                    modification_flag=cur_poly_probe_position_modification_flag,
                    auth_name=cur_poly_probe_position_auth_name,
                    mutated_chem_comp_id=None if (
                            ('Mutated_residue_donor' not in self.list_of_object_indices_refmeas[i].keys())
                            or cur_poly_probe_position_mutation_flag == False)
                    else self.list_mutated_residues[
                        self.list_of_object_indices_refmeas[i]['Mutated_residue_donor']],
                    modified_chem_descriptor=None if (
                                ('Modified_residue_donor' not in self.list_of_object_indices_refmeas[i].keys())
                                or cur_poly_probe_position_modification_flag == False)
                    else self.list_modified_residues[
                        self.list_of_object_indices_refmeas[i]['Modified_residue_donor']])

                cur_poly_probe_position_index = -1
                if not occurs_in_list(cur_poly_probe_position,
                                      self.list_poly_probe_positions):
                    self.list_poly_probe_positions.append(cur_poly_probe_position)
                cur_poly_probe_position_index = self.list_poly_probe_positions.index(cur_poly_probe_position)
                self.list_of_object_indices_refmeas[i]['Poly_probe_position_donor'] = \
                    cur_poly_probe_position_index

        # acceptor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                cur_poly_probe_position_auth_name = None if (
                        'RefMeas_Poly_probe_position_acceptor_name' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_acceptor_name'][i])) else \
                    xls_data['RefMeas_Poly_probe_position_acceptor_name'][i]
                cur_poly_probe_position_entity = \
                    xls_data['RefMeas_Poly_probe_position_acceptor_entity'][i]
                cur_poly_probe_position_asym_unit_id = None if (
                        'RefMeas_Poly_probe_position_acceptor_AsymUnit' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_acceptor_AsymUnit'][i])) else \
                    xls_data['RefMeas_Poly_probe_position_acceptor_AsymUnit'][i]
                cur_poly_probe_position_seq_id = \
                    int(xls_data['RefMeas_Poly_probe_position_acceptor_seq_id'][i])
                cur_poly_probe_position_comp_id = \
                    xls_data['RefMeas_Poly_probe_position_acceptor_comp_id'][i]
                cur_poly_probe_position_atom_id = None if (
                        'RefMeas_Poly_probe_position_acceptor_atom_id' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Poly_probe_position_acceptor_atom_id'][i])) else \
                    xls_data['RefMeas_Poly_probe_position_acceptor_atom_id'][i]
                cur_poly_probe_position_mutation_flag = False if (
                        'RefMeas_Poly_probe_position_acceptor_mutation_flag' not in xls_data.keys()
                        or xls_data['RefMeas_Poly_probe_position_acceptor_mutation_flag'][i] in
                        ['False', 'no', 'No', 'NO']) else True
                cur_poly_probe_position_modification_flag = False if (
                        'RefMeas_Poly_probe_position_acceptor_modification_flag' not in xls_data.keys()
                        or xls_data['RefMeas_Poly_probe_position_acceptor_modification_flag'][i] in
                        ['False', 'no', 'No', 'NO']) else True

                ## Create the Residue or atom
                cur_resatom_new = None
                cur_entity = self.list_ihm_entities[
                    self.list_ihm_entity_ids.index(cur_poly_probe_position_entity)]
                cur_asym_unit = None if cur_poly_probe_position_asym_unit_id is None else \
                    self.list_asym_units[
                        self.list_asym_units_ids.index(cur_poly_probe_position_asym_unit_id)]
                ## First create the residue
                cur_residue = cur_entity.residue(seq_id=cur_poly_probe_position_seq_id)
                if cur_asym_unit is not None:
                    cur_residue.asym = cur_asym_unit

                if cur_poly_probe_position_atom_id is not None:
                    cur_atom = cur_residue.atom(atom_id=cur_poly_probe_position_atom_id)
                    cur_resatom_new = cur_atom
                else:
                    cur_resatom_new = cur_residue
                ## add it to the list of resatoms if it is not there yet
                if get_resatom_from_list(cur_resatom_new, self.list_resatoms) is None:
                    self.list_resatoms.append(cur_resatom_new)
                ## and get the respective entry (to avoid duplicate entries)
                cur_resatom = get_resatom_from_list(cur_resatom_new,
                                                    self.list_resatoms)

                cur_poly_probe_position = ihm.flr.PolyProbePosition(
                    resatom=cur_resatom,
                    mutation_flag=cur_poly_probe_position_mutation_flag,
                    modification_flag=cur_poly_probe_position_modification_flag,
                    auth_name=cur_poly_probe_position_auth_name,
                    mutated_chem_comp_id=None if
                    (('Mutated_residue_acceptor' not in self.list_of_object_indices_refmeas[i].keys())
                     or cur_poly_probe_position_mutation_flag == False) else
                    self.list_mutated_residues[list_of_object_indices_refmeas[i]['Mutated_residue_acceptor']],
                    modified_chem_descriptor=None if (
                                ('Modified_residue_acceptor' not in self.list_of_object_indices_refmeas[i].keys())
                                or cur_poly_probe_position_modification_flag == False) else
                    self.list_modified_residues[list_of_object_indices_refmeas[i]['Modified_residue_acceptor']])

                cur_poly_probe_position_index = -1
                if not occurs_in_list(cur_poly_probe_position,
                                      self.list_poly_probe_positions):
                    self.list_poly_probe_positions.append(cur_poly_probe_position)
                cur_poly_probe_position_index = \
                    self.list_poly_probe_positions.index(cur_poly_probe_position)
                self.list_of_object_indices_refmeas[i]['Poly_probe_position_acceptor'] = cur_poly_probe_position_index

        # => Sample_probe_details
        # donor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i]:
                cur_sample_probe_details_donor_description = None if (
                        'RefMeas_Sample_probe_details_donor_description' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Sample_probe_details_donor_description'][i])) \
                    else xls_data['RefMeas_Sample_probe_details_donor_description'][i]
                cur_sample_probe_details = ihm.flr.SampleProbeDetails(
                    sample=self.list_samples[self.list_of_object_indices_refmeas[i]['Sample']],
                    probe=self.list_probe_donors[self.list_of_object_indices_refmeas[i]['Probe_donor']],
                    fluorophore_type='donor',
                    poly_probe_position=self.list_poly_probe_positions[
                        self.list_of_object_indices_refmeas[i]['Poly_probe_position_donor']],
                    description=cur_sample_probe_details_donor_description)
                cur_sample_probe_details_index = -1
                if cur_sample_probe_details not in self.list_sample_probe_details:
                    self.list_sample_probe_details.append(cur_sample_probe_details)
                cur_sample_probe_details_index = \
                    self.list_sample_probe_details.index(cur_sample_probe_details)
                self.list_of_object_indices_refmeas[i]['Sample_probe_details_donor'] = \
                    cur_sample_probe_details_index

        # acceptor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                cur_sample_probe_details_acceptor_description = None if (
                        'RefMeas_Sample_probe_details_acceptor_description' not in
                        xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Sample_probe_details_acceptor_description'][i])) \
                    else xls_data['RefMeas_Sample_probe_details_acceptor_description'][i]
                cur_sample_probe_details = ihm.flr.SampleProbeDetails(
                    sample=self.list_samples[self.list_of_object_indices_refmeas[i]['Sample']],
                    probe=self.list_probe_acceptors[self.list_of_object_indices_refmeas[i]['Probe_acceptor']],
                    fluorophore_type='acceptor',
                    poly_probe_position=self.list_poly_probe_positions[
                        self.list_of_object_indices_refmeas[i]['Poly_probe_position_acceptor']],
                    description=cur_sample_probe_details_acceptor_description)
                cur_sample_probe_details_index = -1
                if cur_sample_probe_details not in self.list_sample_probe_details:
                    self.list_sample_probe_details.append(cur_sample_probe_details)
                cur_sample_probe_details_index = \
                    self.list_sample_probe_details.index(cur_sample_probe_details)
                self.list_of_object_indices_refmeas[i]['Sample_probe_details_acceptor'] = \
                    cur_sample_probe_details_index

        # => Probe_descriptor_conjugate
        # donor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i]:
                cur_probe_conjugate_descriptor_donor = ihm.ChemDescriptor(
                    auth_name=xls_data['RefMeas_Probe_conjugate_donor_name'][i],
                    chem_comp_id=None if (
                            'RefMeas_Probe_conjugate_donor_chem_comp_id' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_chem_comp_id'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_chem_comp_id'][i],
                    chemical_name=None if (
                            'RefMeas_Probe_conjugate_donor_chemical_name' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_chemical_name'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_chemical_name'][i],
                    common_name=None if (
                            'RefMeas_Probe_conjugate_donor_common_name' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_common_name'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_common_name'][i],
                    smiles=None if (
                            'RefMeas_Probe_conjugate_donor_smiles' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_smiles'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_smiles'][i],
                    smiles_canonical=None if (
                            'RefMeas_Probe_conjugate_donor_smiles_canonical' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_smiles_canonical'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_smiles_canonical'][i],
                    inchi=None if (
                            'RefMeas_Probe_conjugate_donor_inchi' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_inchi'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_inchi'][i],
                    inchi_key=None if (
                            'RefMeas_Probe_conjugate_donor_inchi_key' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_donor_inchi_key'][i]))
                    else xls_data['RefMeas_Probe_conjugate_donor_inchi_key'][i])

                if not occurs_in_list(
                        cur_probe_conjugate_descriptor_donor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(cur_probe_conjugate_descriptor_donor)
                else:
                    cur_probe_conjugate_descriptor_donor = \
                        [x for x in self.list_chemical_descriptors if
                         x.__dict__ == cur_probe_conjugate_descriptor_donor.__dict__][0]

                cur_probe_conjugate_descriptor_donor_index = -1
                if not occurs_in_list(
                        cur_probe_conjugate_descriptor_donor,
                        self.list_probe_conjugate_descriptors):
                    self.list_probe_conjugate_descriptors.append(cur_probe_conjugate_descriptor_donor)
                cur_probe_conjugate_descriptor_donor_index = \
                    self.list_probe_conjugate_descriptors.index(cur_probe_conjugate_descriptor_donor)
                self.list_of_object_indices_refmeas[i]['Probe_conjugate_descriptor_donor'] = \
                    cur_probe_conjugate_descriptor_donor_index

        # acceptor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                cur_probe_conjugate_descriptor_acceptor = ihm.ChemDescriptor(
                    auth_name=xls_data['RefMeas_Probe_conjugate_acceptor_name'][i],
                    chem_comp_id=None if (
                            'RefMeas_Probe_conjugate_acceptor_chem_comp_id' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_acceptor_chem_comp_id'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_chem_comp_id'][i],
                    chemical_name=None if (
                            'RefMeas_Probe_conjugate_acceptor_chemical_name' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_acceptor_chemical_name'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_chemical_name'][i],
                    common_name=None if (
                            'RefMeas_Probe_conjugate_acceptor_common_name' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_acceptor_common_name'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_common_name'][i],
                    smiles=None if (
                            'RefMeas_Probe_conjugate_acceptor_smiles' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_acceptor_smiles'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_smiles'][i],
                    smiles_canonical=None if (
                            'RefMeas_Probe_conjugate_acceptor_smiles_canonical' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_acceptor_smiles_canonical'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_smiles_canonical'][i],
                    inchi=None if (
                            'RefMeas_Probe_conjugate_acceptor_inchi' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['RefMeas_Probe_conjugate_acceptor_inchi'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_inchi'][i],
                    inchi_key=None if (
                            'RefMeas_Probe_conjugate_acceptor_inchi_key' not in
                            xls_data.keys() or pandas.isnull(
                                xls_data['RefMeas_Probe_conjugate_acceptor_inchi_key'][i]))
                    else xls_data['RefMeas_Probe_conjugate_acceptor_inchi_key'][i])

                if not occurs_in_list(
                        cur_probe_conjugate_descriptor_acceptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(cur_probe_conjugate_descriptor_acceptor)
                else:
                    cur_probe_conjugate_descriptor_acceptor = [
                        x for x in self.list_chemical_descriptors
                        if x.__dict__ == cur_probe_conjugate_descriptor_acceptor.__dict__][0]

                cur_probe_conjugate_descriptor_acceptor_index = -1
                if not occurs_in_list(
                        cur_probe_conjugate_descriptor_acceptor,
                        self.list_probe_conjugate_descriptors):
                    self.list_probe_conjugate_descriptors.append(cur_probe_conjugate_descriptor_acceptor)
                cur_probe_conjugate_descriptor_acceptor_index = \
                    self.list_probe_conjugate_descriptors.index(cur_probe_conjugate_descriptor_acceptor)
                self.list_of_object_indices_refmeas[i]['Probe_conjugate_descriptor_acceptor'] = \
                    cur_probe_conjugate_descriptor_acceptor_index

        # => Poly_probe_conjugate
        # donor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_donor_info[i]:
                cur_poly_probe_conjugate_ambiguous_stoichiometry_donor = False if (
                        'RefMeas_Probe_conjugate_ambiguous_probe_stoichiometry_donor' not in xls_data.keys()
                        or xls_data['RefMeas_Probe_conjugate_ambiguous_probe_stoichiometry_donor'][i] in
                        ['False', 'no', 'No', 'NO']) else True
                cur_poly_probe_conjugate_probe_stoichiometry_donor = None if (
                        'RefMeas_Probe_conjugate_probe_stoichiometry_donor' not in xls_data.keys()
                        or pandas.isnull(xls_data['RefMeas_Probe_conjugate_probe_stoichiometry_donor'][i])) \
                    else xls_data['RefMeas_Probe_conjugate_probe_stoichiometry_donor'][i]
                cur_poly_probe_conjugate_donor = ihm.flr.PolyProbeConjugate(
                    sample_probe=self.list_sample_probe_details[
                        self.list_of_object_indices_refmeas[i]['Sample_probe_details_donor']],
                    chem_descriptor=self.list_probe_conjugate_descriptors[
                        self.list_of_object_indices_refmeas[i]['Probe_conjugate_descriptor_donor']],
                    ambiguous_stoichiometry=cur_poly_probe_conjugate_ambiguous_stoichiometry_donor,
                    probe_stoichiometry=cur_poly_probe_conjugate_probe_stoichiometry_donor)
                cur_poly_probe_conjugate_index = -1
                if cur_poly_probe_conjugate_donor not in self.list_poly_probe_conjugates:
                    self.list_poly_probe_conjugates.append(cur_poly_probe_conjugate_donor)
                cur_poly_probe_conjugate_index = \
                    self.list_poly_probe_conjugates.index(cur_poly_probe_conjugate_donor)
                self.list_of_object_indices_refmeas[i]['Poly_probe_conjugate_donor'] = cur_poly_probe_conjugate_index
        # acceptor
        for i in range(nr_entries):
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                cur_poly_probe_conjugate_ambiguous_stoichiometry_acceptor = False if (
                        'RefMeas_Probe_conjugate_ambiguous_probe_stoichiometry_acceptor' not in xls_data.keys()
                        or xls_data['RefMeas_Probe_conjugate_ambiguous_probe_stoichiometry_acceptor'][i]
                        in ['False', 'no', 'No', 'NO']) else True
                cur_poly_probe_conjugate_probe_stoichiometry_acceptor = None if (
                        'RefMeas_Probe_conjugate_probe_stoichiometry_acceptor' not in xls_data.keys()
                        or pandas.isnull(xls_data['RefMeas_Probe_conjugate_probe_stoichiometry_acceptor'][i])) \
                    else xls_data['RefMeas_Probe_conjugate_probe_stoichiometry_acceptor'][i]
                cur_poly_probe_conjugate_acceptor = ihm.flr.PolyProbeConjugate(
                    sample_probe=self.list_sample_probe_details[
                        self.list_of_object_indices_refmeas[i]['Sample_probe_details_acceptor']],
                    chem_descriptor=self.list_probe_conjugate_descriptors[
                        self.list_of_object_indices_refmeas[i]['Probe_conjugate_descriptor_acceptor']],
                    ambiguous_stoichiometry=cur_poly_probe_conjugate_ambiguous_stoichiometry_acceptor,
                    probe_stoichiometry=cur_poly_probe_conjugate_probe_stoichiometry_acceptor)

                cur_poly_probe_conjugate_index = -1
                if cur_poly_probe_conjugate_acceptor not in self.list_poly_probe_conjugates:
                    self.list_poly_probe_conjugates.append(cur_poly_probe_conjugate_acceptor)
                cur_poly_probe_conjugate_index = \
                    self.list_poly_probe_conjugates.index(cur_poly_probe_conjugate_acceptor)
                self.list_of_object_indices_refmeas[i]['Poly_probe_conjugate_acceptor'] = \
                    cur_poly_probe_conjugate_index

        # => ReferenceMeasurement
        tmp_dict_for_reference_measurement_groups = {}
        for i in range(nr_entries):
            cur_ref_measurement_id = xls_data['RefMeas_ID'][i]
            cur_ref_measurement_details = None if (
                    'RefMeas_details' not in xls_data.keys()
                    or pandas.isnull(xls_data['RefMeas_details'][i])) \
                else xls_data['RefMeas_details'][i]
            cur_ref_measurement_group_id = xls_data['RefMeas_Group_ID'][i]

            cur_ref_measurement = None
            # If donor information is given for the reference measurement
            if tmp_list_ref_measurement_contains_donor_info[i]:
                ## Get the respective sample_probe_details entry
                cur_ref_sample_probe = self.list_sample_probe_details[
                    self.list_of_object_indices_refmeas[i]['Sample_probe_details_donor']]

                # Create the RefMeasurement, handle the lifetimes subsequently
                cur_ref_measurement = ihm.flr.RefMeasurement(
                    ref_sample_probe=cur_ref_sample_probe,
                    details=cur_ref_measurement_details)
                # Handle the lifetime lists
                cur_species_fraction_list = str(xls_data['RefMeas_Lifetime_species_fraction_list'][i])
                cur_lifetime_list = str(xls_data['RefMeas_Lifetime_lifetime_list'][i])
                cur_species_names_list = None if (
                            'RefMeas_Lifetime_species_name_list' not in
                            xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Lifetime_species_name_list'][i])) \
                    else str(xls_data['RefMeas_Lifetime_species_name_list'][i])
                # split the semicolon-separated list
                all_cur_species_fractions = cur_species_fraction_list.split(";") \
                    if ";" in cur_species_fraction_list \
                    else [cur_species_fraction_list]
                all_cur_lifetimes = cur_lifetime_list.split(";") \
                    if ";" in cur_lifetime_list \
                    else [cur_lifetime_list]
                all_cur_species_names = [] if cur_species_names_list is None \
                    else (cur_species_names_list.split(";")
                          if ";" in cur_species_names_list
                          else [cur_species_names_list])
                if len(all_cur_species_fractions) != len(all_cur_lifetimes):
                    print(
                        'ERROR: Number of species-fractions (%i) and '
                        'lifetimes (%i) for reference measurement '
                        '(Tab \'Reference_measurements\') differs.'
                        % (len(all_cur_species_fractions), len(all_cur_lifetimes)))
                # Create new RefMeasurementLifetime objects
                for tmpj in range(len(all_cur_species_fractions)):
                    cur_ref_measurement_lifetime = ihm.flr.RefMeasurementLifetime(
                        species_fraction=float(all_cur_species_fractions[tmpj]),
                        lifetime=float(all_cur_lifetimes[tmpj]),
                        species_name=None if (
                                len(all_cur_species_names) != len(all_cur_species_fractions))
                        else all_cur_species_names[tmpj])
                    if cur_ref_measurement_lifetime not in self.list_ref_measurement_lifetimes:
                        self.list_ref_measurement_lifetimes.append(cur_ref_measurement_lifetime)
                    cur_ref_measurement.add_lifetime(cur_ref_measurement_lifetime)

            # If acceptor information is given for the reference measurement
            if tmp_list_ref_measurement_contains_acceptor_info[i]:
                # Get the respective sample_probe_details entry
                cur_ref_sample_probe = self.list_sample_probe_details[
                    self.list_of_object_indices_refmeas[i]['Sample_probe_details_acceptor']]
                # Create the RefMeasurement, handle the lifetimes subsequently
                cur_ref_measurement = ihm.flr.RefMeasurement(
                    ref_sample_probe=cur_ref_sample_probe,
                    details=cur_ref_measurement_details)
                # Handle the lifetime lists
                cur_species_fraction_list = \
                    xls_data['RefMeas_Lifetime_species_fraction_list'][i]
                cur_lifetime_list = xls_data['RefMeas_Lifetime_lifetime_list'][i]
                cur_species_names_list = None if (
                            'RefMeas_Lifetime_species_name_list' not in
                            xls_data.keys() or pandas.isnull(
                    xls_data['RefMeas_Lifetime_species_name_list'][i])) \
                    else xls_data['RefMeas_Lifetime_species_name_list'][i]
                # split the semicolon-separated list
                all_cur_species_fractions = cur_species_fraction_list.split(";")
                all_cur_lifetimes = cur_lifetime_list.split(";")
                all_cur_species_names = [] if cur_species_names_list is None \
                    else cur_species_names_list.split(";")
                if len(all_cur_species_fractions) != len(all_cur_lifetimes):
                    print(
                        'ERROR: Number of species-fractions (%i) and '
                        'lifetimes (%i) for reference measurement '
                        '(Tab \'Reference_measurements\' differs.'
                        % (len(all_cur_species_fractions), len(all_cur_lifetimes)))
                # Create new RefMeasurementLifetime objects
                for tmpj in range(len(all_cur_species_fractions)):
                    cur_ref_measurement_lifetime = ihm.flr.RefMeasurementLifetime(
                        species_fraction=float(all_cur_species_fractions[tmpj]),
                        lifetime=float(all_cur_lifetimes[tmpj]),
                        species_name=None if (
                                len(all_cur_species_names) != len(all_cur_species_fractions))
                        else all_cur_species_names[tmpj])
                    if cur_ref_measurement_lifetime not in self.list_ref_measurement_lifetimes:
                        self.list_ref_measurement_lifetimes.append(cur_ref_measurement_lifetime)
                    cur_ref_measurement.add_lifetime(cur_ref_measurement_lifetime)

            if cur_ref_measurement not in self.list_ref_measurements:
                self.list_ref_measurements.append(cur_ref_measurement)
            cur_ref_measurement_index = self.list_ref_measurements.index(cur_ref_measurement)
            self.list_of_object_indices_refmeas[i]['Reference_measurement'] = cur_ref_measurement_index
            # Store the reference measurement for the reference measurement group
            if cur_ref_measurement_group_id not in tmp_dict_for_reference_measurement_groups.keys():
                tmp_dict_for_reference_measurement_groups[cur_ref_measurement_group_id] = []
            tmp_dict_for_reference_measurement_groups[
                cur_ref_measurement_group_id].append(cur_ref_measurement)

        for i in range(nr_entries):
            cur_ref_measurement_group_id = xls_data['RefMeas_Group_ID'][i]
            cur_ref_measurement_group_details = None if (
                    'RefMeas_Group_Details' not in xls_data.keys() or
                    pandas.isnull(xls_data['RefMeas_Group_Details'][i])) \
                else xls_data['RefMeas_Group_Details'][i]
            cur_ref_measurement_group = \
                ihm.flr.RefMeasurementGroup(details=cur_ref_measurement_group_details)
            ## add the reference measurements for the current reference measurement group
            for entry in tmp_dict_for_reference_measurement_groups[cur_ref_measurement_group_id]:
                cur_ref_measurement_group.add_ref_measurement(entry)
            if cur_ref_measurement_group not in self.list_ref_measurement_groups:
                self.list_ref_measurement_groups.append(cur_ref_measurement_group)
                self.list_ref_measurement_group_ids.append(cur_ref_measurement_group_id)


    def add_flr(self, sheet_name='FLR', skiprows=3, header=0):
        """Read flrCIF information from the excel sheet"""
        if self.verbose:
            print(' ... Processing tab \'FLR\' ...')
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)

        nr_entries = len(xls_data['FLR_ID'])

        # This list of dictionaries contains the indices for the respective
        # objects that they refer to, since entries might refer to the same
        # objects
        self.list_of_object_indices = []
        for i in range(nr_entries):
            self.list_of_object_indices.append({})

        # => Sample conditions
        for i in range(nr_entries):
            # get the respective column for the sample condition
            cur_sample_condition_details = None if (
                    'FLR_Sample_condition' not in xls_data.keys() or
                    pandas.isnull(xls_data['FLR_Sample_condition'][i])) \
                else xls_data['FLR_Sample_condition'][i]
            # create the sample condition object
            cur_sample_condition = ihm.flr.SampleCondition(
                details=cur_sample_condition_details)
            # check whether it is already in the list
            cur_sample_condition_index = -1
            if cur_sample_condition not in self.list_sample_conditions:
                self.list_sample_conditions.append(cur_sample_condition)
            # and store the index of the respective object
            cur_sample_condition_index = \
                self.list_sample_conditions.index(cur_sample_condition)
            self.list_of_object_indices[i]['Sample_condition'] = cur_sample_condition_index

        # => Instrument
        for i in range(nr_entries):
            ## get the respective column for the sample condition
            cur_instrument_details = None if (
                    'FLR_Instrument' not in xls_data.keys() or
                    pandas.isnull(xls_data['FLR_Instrument'][i])) \
                else xls_data['FLR_Instrument'][i]
            # create the object
            cur_instrument = ihm.flr.Instrument(details=cur_instrument_details)
            # check whether it is already in the list
            cur_instrument_index = -1
            if cur_instrument not in self.list_instruments:
                self.list_instruments.append(cur_instrument)
            # and store the index of the respective object
            cur_instrument_index = self.list_instruments.index(cur_instrument)
            self.list_of_object_indices[i]['Instrument'] = cur_instrument_index

        # => Instrument settings
        for i in range(nr_entries):
            cur_inst_setting_details = None if (
                    'FLR_Instrument_setting' not in xls_data.keys() or
                    pandas.isnull(xls_data['FLR_Instrument_setting'][i])) \
                else xls_data['FLR_Instrument_setting'][i]
            cur_inst_setting = ihm.flr.InstSetting(details=cur_inst_setting_details)
            cur_inst_setting_index = -1
            if cur_inst_setting not in self.list_inst_settings:
                self.list_inst_settings.append(cur_inst_setting)
            cur_inst_setting_index = self.list_inst_settings.index(cur_inst_setting)
            self.list_of_object_indices[i]['Instrument_setting'] = cur_inst_setting_index

        # => Experimental condition
        for i in range(nr_entries):
            cur_exp_condition_details = None if (
                    'FLR_Experimental_condition' not in
                    xls_data.keys() or pandas.isnull(
                xls_data['FLR_Experimental_condition'][i])) \
                else xls_data['FLR_Experimental_condition'][i]
            cur_exp_condition = ihm.flr.ExpCondition(details=cur_exp_condition_details)
            cur_exp_condition_index = -1
            if cur_exp_condition not in self.list_exp_conditions:
                self.list_exp_conditions.append(cur_exp_condition)
            cur_exp_condition_index = self.list_exp_conditions.index(cur_exp_condition)
            self.list_of_object_indices[i]['Experimental_condition'] = cur_exp_condition_index

        # => Samples
        for i in range(nr_entries):
            cur_sample_id = xls_data['FLR_Sample']
            cur_sample_num_of_probes = int(xls_data['FLR_Sample_num_of_probes'][i])
            cur_sample_solvent_phase = xls_data['FLR_Sample_solvent_phase'][i]
            cur_sample_description = None if (
                    'FLR_Sample_description' not in xls_data.keys() or
                    pandas.isnull(xls_data['FLR_Sample_description'][i])) \
                else xls_data['FLR_Sample_description'][i]
            cur_sample_details = None if (
                    'FLR_Sample_details' not in xls_data.keys() or
                    pandas.isnull(xls_data['FLR_Sample_details'][i])) \
                else xls_data['FLR_Sample_details'][i]
            cur_entity_assembly = xls_data['FLR_Entity_assembly'][i]
            # create the object
            # The entity assembly is taken from the list of entity assemblies
            # by getting the entry that corresponds to the given
            # entity_assembly_id (cur_entity_assembly)
            cur_sample = ihm.flr.Sample(
                entity_assembly=self.list_entity_assemblies[
                    self.list_entity_assembly_ids.index(cur_entity_assembly)],
                num_of_probes=cur_sample_num_of_probes,
                condition=self.list_sample_conditions[
                    self.list_of_object_indices[i]['Sample_condition']],
                description=cur_sample_description,
                details=cur_sample_details,
                solvent_phase=cur_sample_solvent_phase)
            if not occurs_in_list(cur_sample, self.list_samples):
                self.list_samples.append(cur_sample)
                self.list_sample_ids.append(cur_sample_id)

            cur_sample_index = self.list_samples.index(cur_sample)
            self.list_of_object_indices[i]['Sample'] = cur_sample_index

        # => Experiment
        # in case of the experiment, there is no need to read something
        # from the file. Everything is defined before
        Experiment_1 = ihm.flr.Experiment()
        cur_experiment_index = -1
        for i in range(nr_entries):
            # make sure that every
            # instrument-inst_setting-exp_condition-sample-combination
            # is there only once
            cur_instrument = self.list_instruments[
                self.list_of_object_indices[i]['Instrument']]
            cur_inst_setting = self.list_inst_settings[
                self.list_of_object_indices[i]['Instrument_setting']]
            cur_exp_condition = self.list_exp_conditions[
                self.list_of_object_indices[i]['Experimental_condition']]
            cur_sample = self.list_samples[self.list_of_object_indices[i]['Sample']]
            cur_experiment_details = None if (
                    'FLR_Experiment_Details' not in xls_data.keys() or
                    pandas.isnull(xls_data['FLR_Experiment_Details'][i])) \
                else xls_data['FLR_Experiment_Details'][i]

            if not Experiment_1.contains(
                    instrument=cur_instrument,
                    inst_setting=cur_inst_setting,
                    exp_condition=cur_exp_condition,
                    sample=cur_sample):
                Experiment_1.add_entry(
                    instrument=cur_instrument,
                    inst_setting=cur_inst_setting,
                    exp_condition=cur_exp_condition,
                    sample=cur_sample,
                    details=cur_experiment_details)
        if Experiment_1 not in self.list_experiments:
            self.list_experiments.append(Experiment_1)
        for i in range(nr_entries):
            cur_experiment_index = self.list_experiments.index(Experiment_1)
            self.list_of_object_indices[i]['Experiment'] = cur_experiment_index

        # => Probe
        # Donor
        # Donor - Probe_list
        for i in range(nr_entries):
            cur_probe_donor_name = xls_data['FLR_Probe_donor_name'][i]
            cur_probe_donor_origin = xls_data['FLR_Probe_donor_origin'][i]
            cur_probe_donor_link_type = xls_data['FLR_Probe_donor_link_type'][i]
            cur_reactive_probe_donor_flag = False
            cur_reactive_probe_donor_name = None
            # TODO! MODIFY TO HANDLE EMPTY CELLS HERE!
            if 'FLR_Reactive_probe_donor_name' in xls_data.keys() and not \
                    pandas.isnull(xls_data['FLR_Reactive_probe_donor_name'][i]):
                cur_reactive_probe_donor_flag = True
                cur_reactive_probe_donor_name = xls_data['FLR_Reactive_probe_donor_name'][i]
            cur_probe_list_donor = ihm.flr.ProbeList(
                chromophore_name=cur_probe_donor_name,
                reactive_probe_flag=cur_reactive_probe_donor_flag,
                reactive_probe_name=cur_reactive_probe_donor_name,
                probe_origin=cur_probe_donor_origin,
                probe_link_type=cur_probe_donor_link_type)
            # Probe_descriptor
            cur_probe_descriptor_donor = None
            cur_reactive_probe_donor_chemical_descriptor = None
            ## index of the reactive_probe_chemical_descriptor in the list_chemical_descriptors

            cur_chromophore_donor_chemical_descriptor = None
            # index of the chromophore in the list_chemical_descriptors
            # Reactive probe donor chemical descriptor
            # if there is a name given and at least a smiles, canonical smiles, inchi, or inchi-key
            if ('FLR_Reactive_probe_donor_name' in xls_data.keys()
                    and not pandas.isnull(xls_data['FLR_Reactive_probe_donor_name'][i])
                    and (('FLR_Reactive_probe_donor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['FLR_Reactive_probe_donor_smiles'][i]))
                         or ('FLR_Reactive_probe_donor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Reactive_probe_donor_smiles_canonical'][i]))
                         or ('FLR_Reactive_probe_donor_inchi' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Reactive_probe_donor_inchi'][i]))
                         or ('FLR_Reactive_probe_donor_inchi_key' in xls_flr_data.keys()
                             and not pandas.isnull(xls_data['FLR_Reactive_probe_donor_inchi_key'][i])))):
                cur_reactive_probe_donor_chemical_descriptor = ihm.ChemDescriptor(
                    auth_name=xls_data['FLR_Reactive_probe_donor_name'][i],
                    chem_comp_id=None if (
                            'FLR_Reactive_probe_donor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_donor_chem_comp_id'][i]))
                    else xls_data['FLR_Reactive_probe_donor_chem_comp_id'][i],
                    chemical_name=None if (
                            'FLR_Reactive_probe_donor_chemical_name' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['FLR_Reactive_probe_donor_chemical_name'][i]))
                    else xls_data['FLR_Reactive_probe_donor_chemical_name'][i],
                    common_name=None if (
                            'FLR_Reactive_probe_donor_common_name' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['FLR_Reactive_probe_donor_common_name'][i]))
                    else xls_data['FLR_Reactive_probe_donor_common_name'][i],
                    smiles=None if (
                            'FLR_Reactive_probe_donor_smiles' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['FLR_Reactive_probe_donor_smiles'][i]))
                    else xls_data['FLR_Reactive_probe_donor_smiles'][i],
                    smiles_canonical=None if (
                            'FLR_Reactive_probe_donor_smiles_canonical' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['FLR_Reactive_probe_donor_smiles_canonical'][i]))
                    else xls_data['FLR_Reactive_probe_donor_smiles_canonical'][i],
                    inchi=None if (
                            'FLR_Reactive_probe_donor_inchi' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['FLR_Reactive_probe_donor_inchi'][i]))
                    else xls_data['FLR_Reactive_probe_donor_inchi'][i],
                    inchi_key=None if (
                            'FLR_Reactive_probe_donor_inchi_key' not in
                            xls_data.keys() or pandas.isnull(
                        xls_data['FLR_Reactive_probe_donor_inchi_key'][i]))
                    else xls_data['FLR_Reactive_probe_donor_inchi_key'][i])

            # Same for the chromophore - name and some structural description
            # (smiles or canonical smiles or inchi or inchi key)
            if ('FLR_Chromophore_donor_name' in xls_data.keys()
                    and not pandas.isnull(xls_data['FLR_Chromophore_donor_name'][i])
                    and (('FLR_Chromophore_donor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['FLR_Chromophore_donor_smiles'][i]))
                    or ('FLR_Chromophore_donor_smiles_canonical' in xls_data.keys()
                        and not pandas.isnull(xls_data['FLR_Chromophore_donor_smiles_canonical'][i]))
                    or ('FLR_Chromophore_donor_inchi' in xls_data.keys()
                        and not pandas.isnull(xls_data['FLR_Chromophore_donor_inchi'][i]))
                    or ('FLR_Chromophore_donor_inchi_key' in xls_data.keys()
                        and not pandas.isnull(xls_data['FLR_Chromophore_donor_inchi_key'][i])))):
                cur_chromophore_donor_chemical_descriptor = ihm.ChemDescriptor(
                    auth_name=xls_data['FLR_Chromophore_donor_name'][i],
                    chem_comp_id=None if (
                            'FLR_Chromophore_donor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_chem_comp_id'][i]))
                    else xls_data['FLR_Chromophore_donor_chem_comp_id'][i],
                    chemical_name=None if (
                            'FLR_Chromophore_donor_chemical_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_chemical_name'][i]))
                    else xls_data['FLR_Chromophore_donor_chemical_name'][i],
                    common_name=None if (
                            'FLR_Chromophore_donor_common_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_common_name'][i]))
                    else xls_data['FLR_Chromophore_donor_common_name'][i],
                    smiles=None if (
                            'FLR_Chromophore_donor_smiles' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_smiles'][i]))
                    else xls_data['FLR_Chromophore_donor_smiles'][i],
                    smiles_canonical=None if (
                            'FLR_Chromophore_donor_smiles_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_smiles_canonical'][i]))
                    else xls_data['FLR_Chromophore_donor_smiles_canonical'][i],
                    inchi=None if (
                            'FLR_Chromophore_donor_inchi' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_inchi'][i]))
                    else xls_data['FLR_Chromophore_donor_inchi'][i],
                    inchi_key=None if (
                            'FLR_Chromophore_donor_inchi_key' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_donor_inchi_key'][i]))
                    else xls_data['FLR_Chromophore_donor_inchi_key'][i])

            # In case the current chemical descriptor is already present in
            # the chemical descriptors list, we use that one instead of
            # creating a new one.
            if cur_reactive_probe_donor_chemical_descriptor is not None:
                if not occurs_in_list(
                        cur_reactive_probe_donor_chemical_descriptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(
                        cur_reactive_probe_donor_chemical_descriptor)
                else:
                    cur_reactive_probe_donor_chemical_descriptor = \
                        [x for x in self.list_chemical_descriptors
                         if x.__dict__ == cur_reactive_probe_donor_chemical_descriptor.__dict__][0]
            if cur_chromophore_donor_chemical_descriptor is not None:
                if not occurs_in_list(
                        cur_chromophore_donor_chemical_descriptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(cur_chromophore_donor_chemical_descriptor)
                else:
                    cur_chromophore_donor_chemical_descriptor = \
                        [x for x in self.list_chemical_descriptors
                         if x.__dict__ == cur_chromophore_donor_chemical_descriptor.__dict__][0]

            cur_probe_descriptor_donor = ihm.flr.ProbeDescriptor(
                reactive_probe_chem_descriptor=cur_reactive_probe_donor_chemical_descriptor,
                chromophore_chem_descriptor=cur_chromophore_donor_chemical_descriptor,
                chromophore_center_atom=None if (
                        'FLR_Chromophore_donor_center_atom' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Chromophore_donor_center_atom'][i]))
                else xls_data['FLR_Chromophore_donor_center_atom'][i])

            if not occurs_in_list(
                    cur_probe_descriptor_donor,
                    self.list_chemical_descriptors):
                self.list_chemical_descriptors.append(cur_probe_descriptor_donor)
            else:
                cur_probe_descriptor_donor = \
                [x for x in self.list_chemical_descriptors
                 if x.__dict__ == cur_probe_descriptor_donor.__dict__][0]

            cur_probe_donor = ihm.flr.Probe(
                probe_list_entry=cur_probe_list_donor,
                probe_descriptor=cur_probe_descriptor_donor)

            cur_probe_donor_index = -1
            if cur_probe_donor not in self.list_probe_donors:
                self.list_probe_donors.append(cur_probe_donor)
            cur_probe_donor_index = self.list_probe_donors.index(cur_probe_donor)
            self.list_of_object_indices[i]['Probe_donor'] = cur_probe_donor_index

        # Acceptor
        # Acceptor - Probe_list
        for i in range(nr_entries):
            cur_probe_acceptor_name = xls_data['FLR_Probe_acceptor_name'][i]
            cur_probe_acceptor_origin = xls_data['FLR_Probe_acceptor_origin'][i]
            cur_probe_acceptor_link_type = xls_data['FLR_Probe_acceptor_link_type'][i]
            cur_reactive_probe_acceptor_flag = False
            cur_reactive_probe_acceptor_name = None
            if 'FLR_Reactive_probe_acceptor_name' in xls_data.keys() \
                    and not pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_name'][i]):
                cur_reactive_probe_acceptor_flag = True
                cur_reactive_probe_acceptor_name = \
                    xls_data['FLR_Reactive_probe_acceptor_name'][i]
            cur_probe_list_acceptor = ihm.flr.ProbeList(
                chromophore_name=cur_probe_acceptor_name,
                reactive_probe_flag=cur_reactive_probe_acceptor_flag,
                reactive_probe_name=cur_reactive_probe_acceptor_name,
                probe_origin=cur_probe_acceptor_origin,
                probe_link_type=cur_probe_acceptor_link_type)
            # Probe_descriptor
            cur_probe_descriptor_acceptor = None
            cur_reactive_probe_acceptor_chemical_descriptor = None
            # index of the reactive_probe_chemical_descriptor in the
            # list_chemical_descriptors
            cur_chromophore_acceptor_chemical_descriptor = None
            # Reactive probe acceptor chemical descriptor
            # if there is a name given and at least a smiles,
            # canonical smiles, inchi, or inchi-key
            if ('FLR_Reactive_probe_acceptor_name' in xls_data.keys()
                    and not pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_name'][i])
                    and (('FLR_Reactive_probe_acceptor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_smiles'][i]))
                         or ('FLR_Reactive_probe_acceptor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_smiles_canonical'][i]))
                         or ('FLR_Reactive_probe_acceptor_inchi' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_inchi'][i]))
                         or ('FLR_Reactive_probe_acceptor_inchi_key' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_inchi_key'][i])))):
                cur_reactive_probe_acceptor_chemical_descriptor = ihm.ChemDescriptor(
                    auth_name=xls_data['FLR_Reactive_probe_acceptor_name'][i],
                    chem_comp_id=None if (
                            'FLR_Reactive_probe_acceptor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_chem_comp_id'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_chem_comp_id'][i],
                    chemical_name=None if (
                            'FLR_Reactive_probe_acceptor_chemical_name' not in xls_data.keys() or
                            pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_chemical_name'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_chemical_name'][i],
                    common_name=None if (
                            'FLR_Reactive_probe_acceptor_common_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_common_name'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_common_name'][i],
                    smiles=None if (
                            'FLR_Reactive_probe_acceptor_smiles' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_smiles'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_smiles'][i],
                    smiles_canonical=None if (
                            'FLR_Reactive_probe_acceptor_smiles_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_smiles_canonical'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_smiles_canonical'][i],
                    inchi=None if (
                            'FLR_Reactive_probe_acceptor_inchi' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_inchi'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_inchi'][i],
                    inchi_key=None if (
                            'FLR_Reactive_probe_acceptor_inchi_key' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Reactive_probe_acceptor_inchi_key'][i]))
                    else xls_data['FLR_Reactive_probe_acceptor_inchi_key'][i])

            if ('FLR_Chromophore_acceptor_name' in xls_data.keys()
                    and not pandas.isnull(xls_data['FLR_Chromophore_acceptor_name'][i])
                    and (('FLR_Chromophore_acceptor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['FLR_Chromophore_acceptor_smiles'][i]))
                         or ('FLR_Chromophore_acceptor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Chromophore_acceptor_smiles_canonical'][i]))
                         or ('FLR_Chromophore_acceptor_inchi' in xls_flr_data.keys()
                             and not pandas.isnull(xls_data['FLR_Chromophore_acceptor_inchi'][i]))
                         or ('FLR_Chromophore_acceptor_inchi_key' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Chromophore_acceptor_inchi_key'][i])))):
                cur_chromophore_acceptor_chemical_descriptor = ihm.ChemDescriptor(
                    auth_name=xls_data['FLR_Chromophore_acceptor_name'][i],
                    chem_comp_id=None if (
                            'FLR_Chromophore_acceptor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_acceptor_chem_comp_id'][i]))
                    else xls_data['FLR_Chromophore_acceptor_chem_comp_id'][i],
                    chemical_name=None if (
                            'FLR_Chromophore_acceptor_chemical_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_acceptor_chemical_name'][i]))
                    else xls_data['FLR_Chromophore_acceptor_chemical_name'][i],
                    common_name=None if (
                                'FLR_Chromophore_acceptor_common_name' not in xls_data.keys()
                                or pandas.isnull(xls_data['FLR_Chromophore_acceptor_common_name'][i]))
                    else xls_data['FLR_Chromophore_acceptor_common_name'][i],
                    smiles=None if (
                            'FLR_Chromophore_acceptor_smiles' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_acceptor_smiles'][i]))
                    else xls_data['FLR_Chromophore_acceptor_smiles'][i],
                    smiles_canonical=None if (
                            'FLR_Chromophore_acceptor_smiles_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_acceptor_smiles_canonical'][i]))
                    else xls_data['FLR_Chromophore_acceptor_smiles_canonical'][i],
                    inchi=None if (
                            'FLR_Chromophore_acceptor_inchi' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_acceptor_inchi'][i]))
                    else xls_data['FLR_Chromophore_acceptor_inchi'][i],
                    inchi_key=None if (
                            'FLR_Chromophore_acceptor_inchi_key' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Chromophore_acceptor_inchi_key'][i]))
                    else xls_data['FLR_Chromophore_acceptor_inchi_key'][i])

            if cur_reactive_probe_acceptor_chemical_descriptor is not None:
                if not occurs_in_list(
                        cur_reactive_probe_acceptor_chemical_descriptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(
                        cur_reactive_probe_acceptor_chemical_descriptor)
                else:
                    cur_reactive_probe_acceptor_chemical_descriptor = \
                        [x for x in self.list_chemical_descriptors
                         if x.__dict__ == cur_reactive_probe_acceptor_chemical_descriptor.__dict__][0]
            if cur_chromophore_acceptor_chemical_descriptor is not None:
                if not occurs_in_list(
                        cur_chromophore_acceptor_chemical_descriptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(
                        cur_chromophore_acceptor_chemical_descriptor)
                else:
                    cur_chromophore_acceptor_chemical_descriptor = \
                        [x for x in self.list_chemical_descriptors
                         if x.__dict__ == cur_chromophore_acceptor_chemical_descriptor.__dict__][0]

            cur_probe_descriptor_acceptor = ihm.flr.ProbeDescriptor(
                reactive_probe_chem_descriptor=cur_reactive_probe_acceptor_chemical_descriptor,
                chromophore_chem_descriptor=cur_chromophore_acceptor_chemical_descriptor,
                chromophore_center_atom=None if (
                        'FLR_Chromophore_acceptor_center_atom' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Chromophore_acceptor_center_atom'][i]))
                else xls_data['FLR_Chromophore_acceptor_center_atom'][i])

            if not occurs_in_list(
                    cur_probe_descriptor_acceptor,
                    self.list_chemical_descriptors):
                self.list_chemical_descriptors.append(cur_probe_descriptor_acceptor)
            else:
                cur_probe_descriptor_acceptor = \
                    [x for x in self.list_chemical_descriptors
                     if x.__dict__ == cur_probe_descriptor_acceptor.__dict__][0]

            cur_probe_acceptor = ihm.flr.Probe(
                probe_list_entry=cur_probe_list_acceptor,
                probe_descriptor=cur_probe_descriptor_acceptor)

            cur_probe_acceptor_index = -1
            if cur_probe_acceptor not in self.list_probe_acceptors:
                self.list_probe_acceptors.append(cur_probe_acceptor)
            cur_probe_acceptor_index = self.list_probe_acceptors.index(cur_probe_acceptor)
            self.list_of_object_indices[i]['Probe_acceptor'] = cur_probe_acceptor_index

        # => Modified and mutated residues
        # donor
        for i in range(nr_entries):
            if ('FLR_Modified_residue_donor_name' in xls_data.keys()
                and not pandas.isnull(xls_data['FLR_Modified_residue_donor_name'][i])) \
                    and (('FLR_Modified_residue_donor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['FLR_Modified_residue_donor_smiles'][i]))
                         or ('FLR_Modified_residue_donor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Modified_residue_donor_smiles_canonical'][i]))
                         or ('FLR_Modified_residue_donor_inchi' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Modified_residue_donor_inchi'][i]))
                         or ('FLR_Modified_residue_donor_inchi_key' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Modified_residue_donor_inchi_key'][i]))):
                cur_modified_residue_donor_chemical_descriptor = ihm.ChemDescriptor(
                    auth_name=xls_data['FLR_Modified_residue_donor_name'][i],
                    chem_comp_id=None if (
                            'FLR_Modified_residue_donor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_chem_comp_id'][i]))
                    else xls_data['FLR_Modified_residue_donor_chem_comp_id'][i],
                    chemical_name=None if (
                            'FLR_Modified_residue_donor_chemical_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_chemical_name'][i]))
                    else xls_data['FLR_Modified_residue_donor_chemical_name'][i],
                    common_name=None if (
                            'FLR_Modified_residue_donor_common_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_common_name'][i]))
                    else xls_data['FLR_Modified_residue_donor_common_name'][i],
                    smiles=None if (
                            'FLR_Modified_residue_donor_smiles' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_smiles'][i]))
                    else xls_data['FLR_Modified_residue_donor_smiles'][i],
                    smiles_canonical=None if (
                            'FLR_Modified_residue_donor_smiles_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_smiles_canonical'][i]))
                    else xls_data['FLR_Modified_residue_donor_smiles_canonical'][i],
                    inchi=None if (
                            'FLR_Modified_residue_donor_inchi' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_inchi'][i]))
                    else xls_data['FLR_Modified_residue_donor_inchi'][i],
                    inchi_key=None if (
                            'FLR_Modified_residue_donor_inchi_key' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_donor_inchi_key'][i]))
                    else xls_data['FLR_Modified_residue_donor_inchi_key'][i])
                if not occurs_in_list(
                        cur_modified_residue_donor_chemical_descriptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(
                        cur_modified_residue_donor_chemical_descriptor)
                else:
                    cur_modified_residue_donor_chemical_descriptor = \
                        [x for x in self.list_chemical_descriptors
                         if x.__dict__ == cur_modified_residue_donor_chemical_descriptor.__dict__][0]

                cur_modified_residue_donor_index = -1
                if not occurs_in_list(
                        cur_modified_residue_donor_chemical_descriptor,
                        self.list_modified_residues):
                    self.list_modified_residues.append(
                        cur_modified_residue_donor_chemical_descriptor)
                cur_modified_residue_donor_index = \
                    self.list_modified_residues.index(cur_modified_residue_donor_chemical_descriptor)
                self.list_of_object_indices[i]['Modified_residue_donor'] = \
                    cur_modified_residue_donor_index

        # acceptor
        for i in range(nr_entries):
            if ('FLR_Modified_residue_acceptor_name' in xls_data.keys()
                and not pandas.isnull(xls_data['FLR_Modified_residue_acceptor_name'][i])) \
                    and (('FLR_Modified_residue_acceptor_smiles' in xls_data.keys()
                          and not pandas.isnull(xls_data['FLR_Modified_residue_acceptor_smiles'][i]))
                         or ('FLR_Modified_residue_acceptor_smiles_canonical' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Modified_residue_acceptor_smiles_canonical'][i]))
                         or ('FLR_Modified_residue_acceptor_inchi' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Modified_residue_acceptor_inchi'][i]))
                         or ('FLR_Modified_residue_acceptor_inchi_key' in xls_data.keys()
                             and not pandas.isnull(xls_data['FLR_Modified_residue_acceptor_inchi_key'][i]))):
                cur_modified_residue_acceptor_chemical_descriptor = ihm.ChemDescriptor(
                    auth_name=xls_data['FLR_Modified_residue_acceptor_name'][i],
                    chem_comp_id=None if (
                            'FLR_Modified_residue_acceptor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_chem_comp_id'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_chem_comp_id'][i],
                    chemical_name=None if (
                            'FLR_Modified_residue_acceptor_chemical_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_chemical_name'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_chemical_name'][i],
                    common_name=None if (
                            'FLR_Modified_residue_acceptor_common_name' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_common_name'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_common_name'][i],
                    smiles=None if (
                            'FLR_Modified_residue_acceptor_smiles' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_smiles'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_smiles'][i],
                    smiles_canonical=None if (
                            'FLR_Modified_residue_acceptor_smiles_canonical' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_smiles_canonical'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_smiles_canonical'][i],
                    inchi=None if (
                            'FLR_Modified_residue_acceptor_inchi' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_inchi'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_inchi'][i],
                    inchi_key=None if (
                            'FLR_Modified_residue_acceptor_inchi_key' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Modified_residue_acceptor_inchi_key'][i]))
                    else xls_data['FLR_Modified_residue_acceptor_inchi_key'][i])

                if not occurs_in_list(
                        cur_modified_residue_acceptor_chemical_descriptor,
                        self.list_chemical_descriptors):
                    self.list_chemical_descriptors.append(
                        cur_modified_residue_acceptor_chemical_descriptor)
                else:
                    cur_modified_residue_acceptor_chemical_descriptor = \
                        [x for x in list_chemical_descriptors
                         if x.__dict__ == cur_modified_residue_acceptor_chemical_descriptor.__dict__][0]

                cur_modified_residue_acceptor_index = -1
                if not occurs_in_list(
                        cur_modified_residue_acceptor_chemical_descriptor,
                        self.list_modified_residues):
                    self.list_modified_residues.append(
                        cur_modified_residue_acceptor_chemical_descriptor)
                cur_modified_residue_acceptor_index = self.list_modified_residues.index(
                    cur_modified_residue_acceptor_chemical_descriptor)
                self.list_of_object_indices[i]['Modified_residue_acceptor'] = \
                    cur_modified_residue_acceptor_index

        # => Mutated residues
        # donor
        for i in range(nr_entries):
            if ('FLR_Mutated_residue_donor_ChemComp_ID' in xls_data.keys()
                    and not pandas.isnull(xls_data['FLR_Mutated_residue_donor_ChemComp_ID'][i])):
                cur_mutated_residue_donor_chemcomp_id = xls_data['FLR_Mutated_residue_donor_ChemComp_ID'][i]
                cur_mutated_residue_donor_code = None if (
                        'FLR_Mutated_residue_donor_Code' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Mutated_residue_donor_Code'][i])) \
                    else xls_data['FLR_Mutated_residue_donor_Code'][i]
                cur_mutated_residue_donor_code_canonical = None if (
                        'FLR_Mutated_residue_donor_Code_canonical' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Mutated_residue_donor_Code_canonical'][i])) \
                    else xls_data['FLR_Mutated_residue_donor_Code_canonical'][i]
                cur_mutated_residue_donor_type = None if (
                        'FLR_Mutated_residue_donor_Type' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Mutated_residue_donor_Type'][i])) \
                    else xls_data['FLR_Mutated_residue_donor_Type'][i]
                cur_mutated_residue_donor_chem_comp = None
                if cur_mutated_residue_donor_type == 'Peptide linking':
                    cur_mutated_residue_donor_chem_comp = ihm.PeptideChemComp(
                        id=cur_mutated_residue_donor_chemcomp_id,
                        code=cur_mutated_residue_donor_code,
                        code_canonical=cur_mutated_residue_donor_code_canonical)
                elif cur_mutated_residue_donor_type == 'L-Peptide linking':
                    cur_mutated_residue_donor_chem_comp = ihm.LPeptideChemComp(
                        id=cur_mutated_residue_donor_chemcomp_id,
                        code=cur_mutated_residue_donor_code,
                        code_canonical=cur_mutated_residue_donor_code_canonical)
                elif cur_mutated_residue_donor_type == 'D-Peptide linking':
                    cur_mutated_residue_donor_chem_comp = ihm.DPeptideChemComp(
                        id=cur_mutated_residue_donor_chemcomp_id,
                        code=cur_mutated_residue_donor_code,
                        code_canonical=cur_mutated_residue_donor_code_canonical)
                elif cur_mutated_residue_donor_type == 'DNA linking':
                    cur_mutated_residue_donor_chem_comp = ihm.DNAChemComp(
                        id=cur_mutated_residue_donor_chemcomp_id,
                        code=cur_mutated_residue_donor_code,
                        code_canonical=cur_mutated_residue_donor_code_canonical)
                elif cur_mutated_residue_donor_type == 'RNA linking':
                    cur_mutated_residue_donor_chem_comp = ihm.RNAChemComp(
                        id=cur_mutated_residue_donor_chemcomp_id,
                        code=cur_mutated_residue_donor_code,
                        code_canonical=cur_mutated_residue_donor_code_canonical)
                else:
                    cur_mutated_residue_donor_chem_comp = ihm.ChemComp(
                        id=cur_mutated_residue_donor_chemcomp_id,
                        code=cur_mutated_residue_donor_code,
                        code_canonical=cur_mutated_residue_donor_code_canonical)

                if not occurs_in_list(
                        cur_mutated_residue_donor_chem_comp,
                        self.list_chemical_components):
                    self.list_chemical_components.append(
                        cur_mutated_residue_donor_chem_comp)
                else:
                    cur_mutated_residue_donor_chem_comp = \
                        [x for x in list_chemical_components
                         if x.__dict__ == cur_mutated_residue_donor_chem_comp.__dict__][0]

                cur_mutated_residue_donor_index = -1
                if not occurs_in_list(
                        cur_mutated_residue_donor_chem_comp,
                        self.list_mutated_residues):
                    self.list_mutated_residues.append(
                        cur_mutated_residue_donor_chem_comp)
                cur_mutated_residue_donor_index = \
                    self.list_mutated_residues.index(cur_mutated_residue_donor_chem_comp)
                self.list_of_object_indices[i]['Mutated_residue_donor'] = \
                    cur_mutated_residue_donor_index

        # acceptor
        for i in range(nr_entries):
            if ('FLR_Mutated_residue_acceptor_ChemComp_ID' in xls_data.keys()
                    and not pandas.isnull(xls_data['FLR_Mutated_residue_acceptor_ChemComp_ID'][i])):
                cur_mutated_residue_acceptor_chemcomp_id = xls_data['FLR_Mutated_residue_acceptor_ChemComp_ID'][i]
                cur_mutated_residue_acceptor_code = None if (
                        'FLR_Mutated_residue_acceptor_Code' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Mutated_residue_acceptor_Code'][i])) \
                    else xls_data['FLR_Mutated_residue_acceptor_Code'][i]
                cur_mutated_residue_acceptor_code_canonical = None if (
                        'FLR_Mutated_residue_acceptor_Code_canonical' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Mutated_residue_acceptor_Code_canonical'][i])) \
                    else xls_data['FLR_Mutated_residue_acceptor_Code_canonical'][i]
                cur_mutated_residue_acceptor_type = None if (
                        'FLR_Mutated_residue_acceptor_Type' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Mutated_residue_acceptor_Type'][i])) \
                    else xls_data['FLR_Mutated_residue_acceptor_Type'][i]
                cur_mutated_residue_acceptor_chem_comp = None
                if cur_mutated_residue_acceptor_type == 'Peptide linking':
                    cur_mutated_residue_acceptor_chem_comp = ihm.PeptideChemComp(
                        id=cur_mutated_residue_acceptor_chemcomp_id,
                        code=cur_mutated_residue_acceptor_code,
                        code_canonical=cur_mutated_residue_acceptor_code_canonical)
                elif cur_mutated_residue_acceptor_type == 'L-Peptide linking':
                    cur_mutated_residue_acceptor_chem_comp = ihm.LPeptideChemComp(
                        id=cur_mutated_residue_acceptor_chemcomp_id,
                        code=cur_mutated_residue_acceptor_code,
                        code_canonical=cur_mutated_residue_acceptor_code_canonical)
                elif cur_mutated_residue_acceptor_type == 'D-Peptide linking':
                    cur_mutated_residue_acceptor_chem_comp = ihm.DPeptideChemComp(
                        id=cur_mutated_residue_acceptor_chemcomp_id,
                        code=cur_mutated_residue_acceptor_code,
                        code_canonical=cur_mutated_residue_acceptor_code_canonical)
                elif cur_mutated_residue_acceptor_type == 'DNA linking':
                    cur_mutated_residue_acceptor_chem_comp = ihm.DNAChemComp(
                        id=cur_mutated_residue_acceptor_chemcomp_id,
                        code=cur_mutated_residue_acceptor_code,
                        code_canonical=cur_mutated_residue_acceptor_code_canonical)
                elif cur_mutated_residue_acceptor_type == 'RNA linking':
                    cur_mutated_residue_acceptor_chem_comp = ihm.RNAChemComp(
                        id=cur_mutated_residue_acceptor_chemcomp_id,
                        code=cur_mutated_residue_acceptor_code,
                        code_canonical=cur_mutated_residue_acceptor_code_canonical)
                else:
                    cur_mutated_residue_acceptor_chem_comp = ihm.ChemComp(
                        id=cur_mutated_residue_acceptor_chemcomp_id,
                        code=cur_mutated_residue_acceptor_code,
                        code_canonical=cur_mutated_residue_acceptor_code_canonical)

                if not occurs_in_list(
                        cur_mutated_residue_acceptor_chem_comp,
                        self.list_chemical_components):
                    self.list_chemical_components.append(
                        cur_mutated_residue_acceptor_chem_comp)
                else:
                    cur_mutated_residue_acceptor_chem_comp = \
                        [x for x in self.list_chemical_components
                         if x.__dict__ == cur_mutated_residue_acceptor_chem_comp.__dict__][0]

                cur_mutated_residue_acceptor_index = -1
                if not occurs_in_list(
                        cur_mutated_residue_acceptor_chem_comp,
                        self.list_mutated_residues):
                    self.list_mutated_residues.append(cur_mutated_residue_acceptor_chem_comp)
                cur_mutated_residue_acceptor_index = \
                    self.list_mutated_residues.index(cur_mutated_residue_acceptor_chem_comp)
                self.list_of_object_indices[i]['Mutated_residue_acceptor'] = \
                    cur_mutated_residue_acceptor_index

        # => Poly_probe_position
        # donor
        for i in range(nr_entries):
            cur_poly_probe_position_auth_name = None if (
                    'FLR_Poly_probe_position_donor_name' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Poly_probe_position_donor_name'][i])) \
                else xls_data['FLR_Poly_probe_position_donor_name'][i]
            cur_poly_probe_position_entity = \
                xls_data['FLR_Poly_probe_position_donor_entity'][i]
            cur_poly_probe_position_asym_unit_id = None if (
                    'FLR_Poly_probe_position_donor_AsymUnit' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Poly_probe_position_donor_AsymUnit'][i])) \
                else xls_data['FLR_Poly_probe_position_donor_AsymUnit'][i]
            cur_poly_probe_position_seq_id = \
                int(xls_data['FLR_Poly_probe_position_donor_seq_id'][i])
            cur_poly_probe_position_comp_id = \
                xls_data['FLR_Poly_probe_position_donor_comp_id'][i]
            cur_poly_probe_position_atom_id = None if (
                    'FLR_Poly_probe_position_donor_atom_id' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Poly_probe_position_donor_atom_id'][i])) \
                else xls_data['FLR_Poly_probe_position_donor_atom_id'][i]
            cur_poly_probe_position_mutation_flag = False if (
                    'FLR_Poly_probe_position_donor_mutation_flag' not in xls_data.keys()
                    or xls_data['FLR_Poly_probe_position_donor_mutation_flag'][i] in
                    ['False', 'no', 'No', 'NO']) else True
            cur_poly_probe_position_modification_flag = False if (
                    'FLR_Poly_probe_position_donor_modification_flag' not in xls_data.keys()
                    or xls_data['FLR_Poly_probe_position_donor_modification_flag'][i] in
                    ['False', 'no', 'No', 'NO']) else True

            # Create the Residue or atom
            cur_resatom_new = None
            cur_entity = self.list_ihm_entities[
                self.list_ihm_entity_ids.index(cur_poly_probe_position_entity)]
            cur_asym_unit = None if cur_poly_probe_position_asym_unit_id is None \
                else self.list_asym_units[self.list_asym_units_ids.index(cur_poly_probe_position_asym_unit_id)]
            # First create the residue
            cur_residue = cur_entity.residue(
                seq_id=cur_poly_probe_position_seq_id)
            if cur_asym_unit is not None:
                cur_residue.asym = cur_asym_unit

            if cur_poly_probe_position_atom_id is not None:
                cur_atom = cur_residue.atom(
                    atom_id=cur_poly_probe_position_atom_id)
                cur_resatom_new = cur_atom
            else:
                cur_resatom_new = cur_residue
            # add it to the list of resatoms if it is not there yet
            if get_resatom_from_list(cur_resatom_new,
                                     self.list_resatoms) is None:
                self.list_resatoms.append(cur_resatom_new)
            # and get the respective entry (to avoid duplicate entries)
            cur_resatom = get_resatom_from_list(cur_resatom_new,
                                                self.list_resatoms)

            # TODO: Modify for chem_comp instead of only the name
            cur_poly_probe_position = ihm.flr.PolyProbePosition(
                resatom=cur_resatom,
                mutation_flag=cur_poly_probe_position_mutation_flag,
                modification_flag=cur_poly_probe_position_modification_flag,
                auth_name=cur_poly_probe_position_auth_name,
                mutated_chem_comp_id=
                None if (('Mutated_residue_donor' not in self.list_of_object_indices[i].keys())
                         or cur_poly_probe_position_mutation_flag == False)
                else self.list_mutated_residues[self.list_of_object_indices[i]['Mutated_residue_donor']],
                modified_chem_descriptor=
                None if (('Modified_residue_donor' not in self.list_of_object_indices[i].keys())
                         or cur_poly_probe_position_modification_flag == False)
                else self.list_modified_residues[self.list_of_object_indices[i]['Modified_residue_donor']])

            cur_poly_probe_position_index = -1
            if not occurs_in_list(
                    cur_poly_probe_position, self.list_poly_probe_positions):
                self.list_poly_probe_positions.append(cur_poly_probe_position)
            cur_poly_probe_position_index = \
                self.list_poly_probe_positions.index(cur_poly_probe_position)
            self.list_of_object_indices[i]['Poly_probe_position_donor'] = \
                cur_poly_probe_position_index

        for i in range(nr_entries):
            cur_poly_probe_position_auth_name = None if (
                    'FLR_Poly_probe_position_acceptor_name' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Poly_probe_position_acceptor_name'][i])) \
                else xls_data['FLR_Poly_probe_position_acceptor_name'][i]
            cur_poly_probe_position_entity = \
                xls_data['FLR_Poly_probe_position_acceptor_entity'][i]
            cur_poly_probe_position_asym_unit_id = None if (
                    'FLR_Poly_probe_position_acceptor_AsymUnit' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Poly_probe_position_acceptor_AsymUnit'][i])) \
                else xls_data['FLR_Poly_probe_position_acceptor_AsymUnit'][i]
            cur_poly_probe_position_seq_id = \
                int(xls_data['FLR_Poly_probe_position_acceptor_seq_id'][i])
            cur_poly_probe_position_comp_id = \
                xls_data['FLR_Poly_probe_position_acceptor_comp_id'][i]
            cur_poly_probe_position_atom_id = None if (
                    'FLR_Poly_probe_position_acceptor_atom_id' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Poly_probe_position_acceptor_atom_id'][i])) \
                else xls_data['FLR_Poly_probe_position_acceptor_atom_id'][i]
            cur_poly_probe_position_mutation_flag = False if (
                    'FLR_Poly_probe_position_acceptor_mutation_flag' not in xls_data.keys()
                    or xls_data['FLR_Poly_probe_position_acceptor_mutation_flag'][i] in
                    ['False', 'no', 'No', 'NO']) else True
            cur_poly_probe_position_modification_flag = False if (
                    'FLR_Poly_probe_position_acceptor_modification_flag' not in xls_data.keys()
                    or xls_data['FLR_Poly_probe_position_acceptor_modification_flag'][i] in
                    ['False', 'no', 'No', 'NO']) else True

            # Create the Residue or atom
            cur_resatom_new = None
            cur_entity = self.list_ihm_entities[
                self.list_ihm_entity_ids.index(cur_poly_probe_position_entity)]
            cur_asym_unit = None if cur_poly_probe_position_asym_unit_id is None \
                else self.list_asym_units[
                self.list_asym_units_ids.index(cur_poly_probe_position_asym_unit_id)]
            # First create the residue
            cur_residue = cur_entity.residue(
                seq_id=cur_poly_probe_position_seq_id)
            if cur_asym_unit is not None:
                cur_residue.asym = cur_asym_unit

            if cur_poly_probe_position_atom_id is not None:
                cur_atom = cur_residue.atom(
                    atom_id=cur_poly_probe_position_atom_id)
                cur_resatom_new = cur_atom
            else:
                cur_resatom_new = cur_residue
            # add it to the list of resatoms if it is not there yet
            if get_resatom_from_list(cur_resatom_new,
                                     self.list_resatoms) is None:
                self.list_resatoms.append(cur_resatom_new)
            ## and get the respective entry (to avoid duplicate entries)
            cur_resatom = get_resatom_from_list(cur_resatom_new,
                                                self.list_resatoms)

            cur_poly_probe_position = ihm.flr.PolyProbePosition(
                resatom=cur_resatom,
                mutation_flag=cur_poly_probe_position_mutation_flag,
                modification_flag=cur_poly_probe_position_modification_flag,
                auth_name=cur_poly_probe_position_auth_name,
                mutated_chem_comp_id=
                None if (('Mutated_residue_acceptor' not in self.list_of_object_indices[i].keys())
                         or cur_poly_probe_position_mutation_flag == False)
                else self.list_mutated_residues[self.list_of_object_indices[i]['Mutated_residue_acceptor']],
                modified_chem_descriptor=
                None if (('Modified_residue_acceptor' not in self.list_of_object_indices[i].keys())
                         or cur_poly_probe_position_modification_flag == False)
                else self.list_modified_residues[self.list_of_object_indices[i]['Modified_residue_acceptor']])

            cur_poly_probe_position_index = -1
            if not occurs_in_list(cur_poly_probe_position,
                                  self.list_poly_probe_positions):
                self.list_poly_probe_positions.append(cur_poly_probe_position)
            cur_poly_probe_position_index = self.list_poly_probe_positions.index(cur_poly_probe_position)
            self.list_of_object_indices[i]['Poly_probe_position_acceptor'] = \
                cur_poly_probe_position_index

        # => Sample_probe_details
        # donor
        for i in range(nr_entries):
            cur_sample_probe_details_donor_description = None if (
                    'FLR_Sample_probe_details_donor_description' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Sample_probe_details_donor_description'][i])) \
                else xls_data['FLR_Sample_probe_details_donor_description'][i]
            cur_sample_probe_details = ihm.flr.SampleProbeDetails(
                sample=self.list_samples[self.list_of_object_indices[i]['Sample']],
                probe=self.list_probe_donors[self.list_of_object_indices[i]['Probe_donor']],
                fluorophore_type='donor',
                poly_probe_position=self.list_poly_probe_positions[
                    self.list_of_object_indices[i]['Poly_probe_position_donor']],
                description=cur_sample_probe_details_donor_description)
            cur_sample_probe_details_index = -1
            if cur_sample_probe_details not in self.list_sample_probe_details:
                self.list_sample_probe_details.append(cur_sample_probe_details)
            cur_sample_probe_details_index = \
                self.list_sample_probe_details.index(cur_sample_probe_details)
            self.list_of_object_indices[i]['Sample_probe_details_donor'] = \
                cur_sample_probe_details_index

        # acceptor
        for i in range(nr_entries):
            cur_sample_probe_details_acceptor_description = None if (
                    'FLR_Sample_probe_details_acceptor_description' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Sample_probe_details_acceptor_description'][i])) \
                else xls_data['FLR_Sample_probe_details_acceptor_description'][i]
            cur_sample_probe_details = ihm.flr.SampleProbeDetails(
                sample=self.list_samples[self.list_of_object_indices[i]['Sample']],
                probe=self.list_probe_acceptors[
                    self.list_of_object_indices[i]['Probe_acceptor']],
                fluorophore_type='acceptor',
                poly_probe_position=self.list_poly_probe_positions[
                    self.list_of_object_indices[i]['Poly_probe_position_acceptor']],
                description=cur_sample_probe_details_acceptor_description)
            cur_sample_probe_details_index = -1
            if cur_sample_probe_details not in self.list_sample_probe_details:
                self.list_sample_probe_details.append(cur_sample_probe_details)
            cur_sample_probe_details_index = \
                self.list_sample_probe_details.index(cur_sample_probe_details)
            self.list_of_object_indices[i]['Sample_probe_details_acceptor'] = \
                cur_sample_probe_details_index

        # => Probe_descriptor_conjugate
        # donor
        for i in range(nr_entries):
            cur_probe_conjugate_descriptor_donor = ihm.ChemDescriptor(
                auth_name=xls_data['FLR_Probe_conjugate_donor_name'][i],
                chem_comp_id=None if (
                            'FLR_Probe_conjugate_donor_chem_comp_id' not in xls_data.keys()
                            or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_chem_comp_id'][i]))
                else xls_data['FLR_Probe_conjugate_donor_chem_comp_id'][i],
                chemical_name=None if (
                        'FLR_Probe_conjugate_donor_chemical_name' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_chemical_name'][i]))
                else xls_data['FLR_Probe_conjugate_donor_chemical_name'][i],
                common_name=None if (
                        'FLR_Probe_conjugate_donor_common_name' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_common_name'][i]))
                else xls_data['FLR_Probe_conjugate_donor_common_name'][i],
                smiles=None if (
                        'FLR_Probe_conjugate_donor_smiles' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_smiles'][i]))
                else xls_data['FLR_Probe_conjugate_donor_smiles'][i],
                smiles_canonical=None if (
                        'FLR_Probe_conjugate_donor_smiles_canonical' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_smiles_canonical'][i]))
                else xls_data['FLR_Probe_conjugate_donor_smiles_canonical'][i],
                inchi=None if (
                        'FLR_Probe_conjugate_donor_inchi' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_inchi'][i]))
                else xls_data['FLR_Probe_conjugate_donor_inchi'][i],
                inchi_key=None if (
                        'FLR_Probe_conjugate_donor_inchi_key' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_donor_inchi_key'][i]))
                else xls_data['FLR_Probe_conjugate_donor_inchi_key'][i])

            if not occurs_in_list(cur_probe_conjugate_descriptor_donor,
                                  self.list_chemical_descriptors):
                self.list_chemical_descriptors.append(cur_probe_conjugate_descriptor_donor)
            else:
                cur_probe_conjugate_descriptor_donor = \
                    [x for x in self.list_chemical_descriptors
                     if x.__dict__ == cur_probe_conjugate_descriptor_donor.__dict__][0]

            cur_probe_conjugate_descriptor_donor_index = -1
            if not occurs_in_list(cur_probe_conjugate_descriptor_donor,
                                  self.list_probe_conjugate_descriptors):
                self.list_probe_conjugate_descriptors.append(cur_probe_conjugate_descriptor_donor)
            cur_probe_conjugate_descriptor_donor_index = \
                self.list_probe_conjugate_descriptors.index(cur_probe_conjugate_descriptor_donor)
            self.list_of_object_indices[i]['Probe_conjugate_descriptor_donor'] = \
                cur_probe_conjugate_descriptor_donor_index

        # acceptor
        for i in range(nr_entries):
            cur_probe_conjugate_descriptor_acceptor = ihm.ChemDescriptor(
                auth_name=xls_data['FLR_Probe_conjugate_acceptor_name'][i],
                chem_comp_id=None if (
                        'FLR_Probe_conjugate_acceptor_chem_comp_id' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_chem_comp_id'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_chem_comp_id'][i],
                chemical_name=None if (
                        'FLR_Probe_conjugate_acceptor_chemical_name' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_chemical_name'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_chemical_name'][i],
                common_name=None if (
                        'FLR_Probe_conjugate_acceptor_common_name' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_common_name'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_common_name'][i],
                smiles=None if (
                        'FLR_Probe_conjugate_acceptor_smiles' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_smiles'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_smiles'][i],
                smiles_canonical=None if (
                        'FLR_Probe_conjugate_acceptor_smiles_canonical' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_smiles_canonical'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_smiles_canonical'][i],
                inchi=None if ('FLR_Probe_conjugate_acceptor_inchi' not in xls_data.keys()
                               or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_inchi'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_inchi'][i],
                inchi_key=None if (
                        'FLR_Probe_conjugate_acceptor_inchi_key' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_Probe_conjugate_acceptor_inchi_key'][i]))
                else xls_data['FLR_Probe_conjugate_acceptor_inchi_key'][i])
            if not occurs_in_list(cur_probe_conjugate_descriptor_acceptor,
                                  self.list_chemical_descriptors):
                self.list_chemical_descriptors.append(
                    cur_probe_conjugate_descriptor_acceptor)
            else:
                cur_probe_conjugate_descriptor_acceptor = \
                    [x for x in self.list_chemical_descriptors
                     if x.__dict__ == cur_probe_conjugate_descriptor_acceptor.__dict__][0]

            cur_probe_conjugate_descriptor_acceptor_index = -1
            if not occurs_in_list(cur_probe_conjugate_descriptor_acceptor,
                                  self.list_probe_conjugate_descriptors):
                self.list_probe_conjugate_descriptors.append(
                    cur_probe_conjugate_descriptor_acceptor)
            cur_probe_conjugate_descriptor_acceptor_index = \
                self.list_probe_conjugate_descriptors.index(cur_probe_conjugate_descriptor_acceptor)
            self.list_of_object_indices[i]['Probe_conjugate_descriptor_acceptor'] = \
                cur_probe_conjugate_descriptor_acceptor_index

        # => Poly_probe_conjugate
        # donor
        for i in range(nr_entries):
            cur_poly_probe_conjugate_ambiguous_stoichiometry_donor = False if (
                    'FLR_Probe_conjugate_ambiguous_probe_stoichiometry_donor' not in xls_data.keys()
                    or xls_data['FLR_Probe_conjugate_ambiguous_probe_stoichiometry_donor'][i] in
                    ['False', 'no', 'No', 'NO']) else True
            cur_poly_probe_conjugate_probe_stoichiometry_donor = None if (
                    'FLR_Probe_conjugate_probe_stoichiometry_donor' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Probe_conjugate_probe_stoichiometry_donor'][i])) \
                else xls_data['FLR_Probe_conjugate_probe_stoichiometry_donor'][i]
            cur_poly_probe_conjugate_donor = ihm.flr.PolyProbeConjugate(
                sample_probe=self.list_sample_probe_details[
                    self.list_of_object_indices[i]['Sample_probe_details_donor']],
                chem_descriptor=self.list_probe_conjugate_descriptors[
                    self.list_of_object_indices[i]['Probe_conjugate_descriptor_donor']],
                ambiguous_stoichiometry=cur_poly_probe_conjugate_ambiguous_stoichiometry_donor,
                probe_stoichiometry=cur_poly_probe_conjugate_probe_stoichiometry_donor)
            cur_poly_probe_conjugate_index = -1
            if cur_poly_probe_conjugate_donor not in self.list_poly_probe_conjugates:
                self.list_poly_probe_conjugates.append(cur_poly_probe_conjugate_donor)
            cur_poly_probe_conjugate_index = \
                self.list_poly_probe_conjugates.index(cur_poly_probe_conjugate_donor)
            self.list_of_object_indices[i]['Poly_probe_conjugate_donor'] = \
                cur_poly_probe_conjugate_index
        # acceptor
        for i in range(nr_entries):
            cur_poly_probe_conjugate_ambiguous_stoichiometry_acceptor = False if (
                    'FLR_Probe_conjugate_ambiguous_probe_stoichiometry_acceptor' not in xls_data.keys() or
                    xls_data['FLR_Probe_conjugate_ambiguous_probe_stoichiometry_acceptor'][i] in
                    ['False', 'no', 'No', 'NO']) else True
            cur_poly_probe_conjugate_probe_stoichiometry_acceptor = None if (
                    'FLR_Probe_conjugate_probe_stoichiometry_acceptor' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Probe_conjugate_probe_stoichiometry_acceptor'][i])) \
                else xls_data['FLR_Probe_conjugate_probe_stoichiometry_acceptor'][i]
            cur_poly_probe_conjugate_acceptor = ihm.flr.PolyProbeConjugate(
                sample_probe=self.list_sample_probe_details[
                    self.list_of_object_indices[i]['Sample_probe_details_acceptor']],
                chem_descriptor=self.list_probe_conjugate_descriptors[
                    self.list_of_object_indices[i]['Probe_conjugate_descriptor_acceptor']],
                ambiguous_stoichiometry=cur_poly_probe_conjugate_ambiguous_stoichiometry_acceptor,
                probe_stoichiometry=cur_poly_probe_conjugate_probe_stoichiometry_acceptor)

            cur_poly_probe_conjugate_index = -1
            if cur_poly_probe_conjugate_acceptor not in self.list_poly_probe_conjugates:
                self.list_poly_probe_conjugates.append(cur_poly_probe_conjugate_acceptor)
            cur_poly_probe_conjugate_index = \
                self.list_poly_probe_conjugates.index(cur_poly_probe_conjugate_acceptor)
            self.list_of_object_indices[i]['Poly_probe_conjugate_acceptor'] = \
                cur_poly_probe_conjugate_index

        # => FRET_Forster_radius
        for i in range(nr_entries):
            cur_forster_radius = float(xls_data['FLR_Forster_radius'][i])
            cur_reduced_forster_radius = float(xls_data['FLR_reduced_Forster_radius'][i])
            cur_forster_radius_object = ihm.flr.FRETForsterRadius(
                donor_probe=self.list_probe_donors[
                    self.list_of_object_indices[i]['Probe_donor']],
                acceptor_probe=self.list_probe_acceptors[
                    self.list_of_object_indices[i]['Probe_acceptor']],
                forster_radius=cur_forster_radius,
                reduced_forster_radius=cur_reduced_forster_radius)

            cur_fret_forster_radius_object_index = -1
            if cur_forster_radius_object not in self.list_fret_forster_radius:
                self.list_fret_forster_radius.append(cur_forster_radius_object)
            cur_fret_forster_radius_object_index = \
                self.list_fret_forster_radius.index(cur_forster_radius_object)
            self.list_of_object_indices[i]['FRET_forster_radius'] = \
                cur_fret_forster_radius_object_index

        # => FRET_calibration_parameters

        for i in range(nr_entries):
            cur_phi_acceptor = None if (
                    'FLR_Calibration_parameters_phi_acceptor' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_phi_acceptor'][i])) \
                else xls_data['FLR_Calibration_parameters_phi_acceptor'][i]
            cur_alpha = None if (
                    'FLR_Calibration_parameters_alpha' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_alpha'][i])) \
                else xls_data['FLR_Calibration_parameters_alpha'][i]
            cur_alpha_sd = None if (
                    'FLR_Calibration_parameters_alpha_sd' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_alpha_sd'][i])) \
                else xls_data['FLR_Calibration_parameters_alpha_sd'][i]
            cur_beta = None if (
                    'FLR_Calibration_parameters_beta' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_beta'][i])) \
                else xls_data['FLR_Calibration_parameters_beta'][i]
            cur_gamma = None if (
                    'FLR_Calibration_parameters_gamma' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_gamma'][i])) \
                else xls_data['FLR_Calibration_parameters_gamma'][i]
            cur_delta = None if (
                    'FLR_Calibration_parameters_delta' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_delta'][i])) \
                else xls_data['FLR_Calibration_parameters_delta'][i]
            cur_gG_gR_ratio = None if (
                    'FLR_Calibration_parameters_gG_gR_ratio' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_gG_gR_ratio'][i])) \
                else xls_data['FLR_Calibration_parameters_gG_gR_ratio'][i]
            cur_a_b = None if (
                    'FLR_Calibration_parameters_a_b' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Calibration_parameters_a_b'][i])) \
                else xls_data['FLR_Calibration_parameters_a_b'][i]
            # Make sure one fret_calibration_parameter is not None
            cur_fret_calibration_parameters = None
            if (cur_phi_acceptor is not None) or (cur_alpha is not None) \
                    or (cur_beta is not None) or (cur_gamma is not None) \
                    or (cur_delta is not None) or (cur_gG_gR_ratio is not None) \
                    or (cur_a_b is not None):
                cur_fret_calibration_parameters = ihm.flr.FRETCalibrationParameters(
                    phi_acceptor=cur_phi_acceptor,
                    alpha=cur_alpha,
                    alpha_sd=cur_alpha_sd,
                    beta=cur_beta,gamma=cur_gamma,
                    delta=cur_delta,
                    gg_gr_ratio=cur_gG_gR_ratio,
                    a_b=cur_a_b)
            else:
                cur_fret_calibration_parameters = None
            cur_fret_calibration_parameters_index = -1
            if cur_fret_calibration_parameters is not None:
                if not occurs_in_list(cur_fret_calibration_parameters,
                                      self.list_fret_calibration_parameters):
                    self.list_fret_calibration_parameters.append(
                        cur_fret_calibration_parameters)
                cur_fret_calibration_parameters_index = \
                    self.list_fret_calibration_parameters.index(cur_fret_calibration_parameters)
                self.list_of_object_indices[i]['FRET_calibration_parameters'] = \
                    cur_fret_calibration_parameters_index

        # => Lifetime Fit Model
        for i in range(nr_entries):
            cur_lifetime_fit_model_name = None if (
                    'FLR_Lifetime_fit_model_name' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Lifetime_fit_model_name'][i])) \
                else xls_data['FLR_Lifetime_fit_model_name'][i]
            cur_lifetime_fit_model_description = None if (
                    'FLR_Lifetime_fit_model_description' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Lifetime_fit_model_description'][i])) \
                else xls_data['FLR_Lifetime_fit_model_description'][i]
            cur_lifetime_fit_model_external_file_id = None if (
                    'FLR_Lifetime_fit_model_external_file' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Lifetime_fit_model_external_file'][i])) \
                else xls_data['FLR_Lifetime_fit_model_external_file'][i]
            cur_lifetime_fit_model_citation_id = None if (
                    'FLR_Lifetime_fit_model_citation' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Lifetime_fit_model_citation'][i])) \
                else xls_data['FLR_Lifetime_fit_model_citation'][i]
            cur_lifetime_fit_model_external_file = None if (
                    cur_lifetime_fit_model_external_file_id is None) \
                else self.list_external_files_locations[
                self.list_external_files_ids.index(cur_lifetime_fit_model_external_file_id)]
            cur_lifetime_fit_model_citation = None if (
                    cur_lifetime_fit_model_citation_id is None) \
                else self.list_citations[
                self.list_citation_ids.index(cur_lifetime_fit_model_citation_id)]
            cur_lifetime_fit_model = ihm.flr.LifetimeFitModel(
                name=cur_lifetime_fit_model_name,
                description=cur_lifetime_fit_model_description,
                external_file=cur_lifetime_fit_model_external_file,
                citation=cur_lifetime_fit_model_citation)
            if cur_lifetime_fit_model not in self.list_lifetime_fit_models:
                self.list_lifetime_fit_models.append(cur_lifetime_fit_model)
            cur_lifetime_fit_model_index = \
                self.list_lifetime_fit_models.index(cur_lifetime_fit_model)
            self.list_of_object_indices[i]['Lifetime_fit_model'] = \
                cur_lifetime_fit_model_index

        # => FRET_analysis
        for i in range(nr_entries):
            cur_fret_analysis_type = None if (
                    'FLR_FRET_analysis_type' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_analysis_type'][i])) \
                else xls_data['FLR_FRET_analysis_type'][i]
            cur_fret_analysis_method_name = None if (
                    'FLR_FRET_analysis_method_name' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_analysis_method_name'][i])) \
                else xls_data['FLR_FRET_analysis_method_name'][i]
            cur_fret_analysis_chi_square_reduced = None if (
                    'FLR_FRET_analysis_chi_squared_reduced' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_analysis_chi_squared_reduced'][i])) \
                else float(xls_data['FLR_FRET_analysis_chi_squared_reduced'][i])
            cur_fret_analysis_experiment_id = self.list_experiments[
                self.list_of_object_indices[i]['Experiment']]
            cur_fret_analysis_sample_probe_id_1 = self.list_sample_probe_details[
                self.list_of_object_indices[i]['Sample_probe_details_donor']]
            cur_fret_analysis_sample_probe_id_2 = self.list_sample_probe_details[
                self.list_of_object_indices[i]['Sample_probe_details_acceptor']]
            cur_fret_analysis_forster_radius_id = self.list_fret_forster_radius[
                self.list_of_object_indices[i]['FRET_forster_radius']]
            cur_fret_analysis_calibration_parameters_id = None if (
                    cur_fret_analysis_type != 'intensity-based') \
                else self.list_fret_calibration_parameters[
                self.list_of_object_indices[i]['FRET_calibration_parameters']]
            cur_lifetime_fit_model_id = None if (
                    'FLR_Lifetime_fit_model_name' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Lifetime_fit_model_name'][i])) \
                else xls_data['FLR_Lifetime_fit_model_name'][i]
            cur_ref_measurement_group_id = None if (
                    'FLR_Ref_measurement_group_id' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Ref_measurement_group_id'][i])) \
                else xls_data['FLR_Ref_measurement_group_id'][i]
            cur_fret_analysis_dataset_list_id = \
                xls_data['FLR_FRET_analysis_dataset_list_id'][i]
            cur_fret_analysis_external_file_id = None if (
                    'FLR_FRET_analysis_external_file_id' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_analysis_external_file_id'][i])) \
                else xls_data['FLR_FRET_analysis_external_file_id'][i]
            cur_fret_analysis_software_id = None if (
                    'FLR_FRET_analysis_software_id' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_analysis_software_id'][i])) \
                else xls_data['FLR_FRET_analysis_software_id'][i]
            cur_lifetime_fit_model = None if (cur_lifetime_fit_model_id is None) \
                else self.list_lifetime_fit_models[
                self.list_of_object_indices[i]['Lifetime_fit_model']]
            cur_ref_measurement_group = None if (cur_ref_measurement_group_id is None) \
                else self.list_ref_measurement_groups[
                self.list_ref_measurement_group_ids.index(cur_ref_measurement_group_id)]
            # get the dataset, external file, and software
            cur_fret_analysis_dataset = self.list_datasets[
                self.list_dataset_ids.index(cur_fret_analysis_dataset_list_id)]
            cur_fret_analysis_external_file = None if (cur_fret_analysis_external_file_id == None) \
                else self.list_external_files_locations[
                self.list_external_files_ids.index(cur_fret_analysis_external_file_id)]
            cur_fret_analysis_software = None if (cur_fret_analysis_software_id == None) \
                else self.list_ihm_softwares[
                self.list_ihm_software_ids.index(cur_fret_analysis_software_id)]
            cur_fret_analysis = ihm.flr.FRETAnalysis(
                experiment=cur_fret_analysis_experiment_id,
                sample_probe_1=cur_fret_analysis_sample_probe_id_1,
                sample_probe_2=cur_fret_analysis_sample_probe_id_2,
                forster_radius=cur_fret_analysis_forster_radius_id,
                type=cur_fret_analysis_type,
                calibration_parameters=cur_fret_analysis_calibration_parameters_id,
                lifetime_fit_model=cur_lifetime_fit_model,
                ref_measurement_group=cur_ref_measurement_group,
                method_name=cur_fret_analysis_method_name,
                chi_square_reduced=cur_fret_analysis_chi_square_reduced,
                dataset=cur_fret_analysis_dataset,
                external_file=cur_fret_analysis_external_file,
                software=cur_fret_analysis_software)

            cur_fret_analysis_index = -1
            if cur_fret_analysis not in self.list_fret_analysis:
                self.list_fret_analysis.append(cur_fret_analysis)
            cur_fret_analysis_index = self.list_fret_analysis.index(cur_fret_analysis)
            self.list_of_object_indices[i]['FRET_analysis'] = cur_fret_analysis_index

        # => Peak_assignment
        for i in range(nr_entries):
            cur_peak_assignment_method_name = None if (
                    'FLR_Peak_assignment_method_name' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Peak_assignment_method_name'][i])) \
                else xls_data['FLR_Peak_assignment_method_name'][i]
            cur_peak_assignment_details = None if (
                    'FLR_Peak_assignment_details' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_Peak_assignment_details'][i])) \
                else xls_data['FLR_Peak_assignment_details'][i]
            cur_peak_assignment = None if (cur_peak_assignment_method_name is None) \
                else ihm.flr.PeakAssignment(
                method_name=cur_peak_assignment_method_name,
                details=cur_peak_assignment_details)
            cur_peak_assignment_index = -1
            if cur_peak_assignment is not None:
                if cur_peak_assignment not in self.list_peak_assignments:
                    self.list_peak_assignments.append(cur_peak_assignment)
                cur_peak_assignment_index = \
                    self.list_peak_assignments.index(cur_peak_assignment)
            self.list_of_object_indices[i]['Peak_assignment'] = \
                cur_peak_assignment_index

        ###### Fret_distance_restraint_group
        for i in range(nr_entries):
            cur_fret_distance_restraint_group_id = None if (
                    'FLR_FRET_distance_restraint_group' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_distance_restraint_group'][i])) \
                else xls_data['FLR_FRET_distance_restraint_group'][i]
            cur_fret_distance_restraint_group = ihm.flr.FRETDistanceRestraintGroup()
            cur_fret_distance_restraint_group._tmp_id = cur_fret_distance_restraint_group_id
            cur_fret_distance_restraint_group_index = -1
            if cur_fret_distance_restraint_group not in self.list_fret_distance_restraint_groups:
                self.list_fret_distance_restraint_groups.append(cur_fret_distance_restraint_group)
            cur_fret_distance_restraint_group_index = \
                self.list_fret_distance_restraint_groups.index(cur_fret_distance_restraint_group)
            self.list_of_object_indices[i]['FRET_distance_restraint_group'] = \
                cur_fret_distance_restraint_group_index

        # => FRET_distance_restraint
        for i in range(nr_entries):
            cur_fret_distance_restraint_id = xls_data['FLR_FRET_distance_restraint_id'][i]
            cur_fret_distance_restraint_state_id = None if (
                    'FLR_FRET_distance_restraint_state_id' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_distance_restraint_state_id'][i])) \
                else self.list_models_states[self.list_models_state_ids.index(
                xls_data['FLR_FRET_distance_restraint_state_id'][i])]
            cur_fret_distance_restraint_distance = \
                float(xls_data['FLR_FRET_distance_restraint_distance'][i])
            cur_fret_distance_restraint_distance_error_plus = None if (
                    'FLR_FRET_distance_restraint_distance_error_plus' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_distance_restraint_distance_error_plus'][i])) \
                else xls_data['FLR_FRET_distance_restraint_distance_error_plus'][i]
            cur_fret_distance_restraint_distance_error_minus = None if (
                    'FLR_FRET_distance_restraint_distance_error_minus' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_distance_restraint_distance_error_minus'][i])) \
                else float(xls_data['FLR_FRET_distance_restraint_distance_error_minus'][i])
            cur_fret_distance_restraint_distance_type = \
                xls_data['FLR_FRET_distance_restraint_distance_type'][i]
            cur_fret_distance_restraint_population_fraction = None if (
                    'FLR_FRET_distance_restraint_population_fraction' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_distance_restraint_population_fraction'][i])) \
                else float(xls_data['FLR_FRET_distance_restraint_population_fraction'][i])
            cur_fret_distance_restraint_sample_probe_id_1 = self.list_sample_probe_details[
                self.list_of_object_indices[i]['Sample_probe_details_donor']]
            cur_fret_distance_restraint_sample_probe_id_2 = self.list_sample_probe_details[
                self.list_of_object_indices[i]['Sample_probe_details_acceptor']]
            cur_fret_distance_restraint_peak_assignment_id = None if (
                    self.list_of_object_indices[i]['Peak_assignment'] == -1) \
                else self.list_peak_assignments[self.list_of_object_indices[i]['Peak_assignment']]
            cur_fret_distance_restraint_analysis_id = self.list_fret_analysis[
                self.list_of_object_indices[i]['FRET_analysis']]
            cur_fret_distance_restraint = ihm.flr.FRETDistanceRestraint(
                sample_probe_1=cur_fret_distance_restraint_sample_probe_id_1,
                sample_probe_2=cur_fret_distance_restraint_sample_probe_id_2,
                state=cur_fret_distance_restraint_state_id,
                analysis=cur_fret_distance_restraint_analysis_id,
                distance=cur_fret_distance_restraint_distance,
                distance_error_plus=cur_fret_distance_restraint_distance_error_plus,
                distance_error_minus=cur_fret_distance_restraint_distance_error_minus,
                distance_type=cur_fret_distance_restraint_distance_type,
                population_fraction=cur_fret_distance_restraint_population_fraction,
                peak_assignment=cur_fret_distance_restraint_peak_assignment_id)

            cur_fret_distance_restraint_index = -1
            if cur_fret_distance_restraint not in self.list_fret_distance_restraints:
                self.list_fret_distance_restraints.append(cur_fret_distance_restraint)
                self.list_fret_distance_restraint_ids.append(cur_fret_distance_restraint_id)
                # add the distance restraint to the distance restraint group
                self.list_fret_distance_restraint_groups[
                    self.list_of_object_indices[i][
                        'FRET_distance_restraint_group']].add_distance_restraint(
                    cur_fret_distance_restraint)
            cur_fret_distance_restraint_index = \
                self.list_fret_distance_restraints.index(cur_fret_distance_restraint)
            self.list_of_object_indices[i]['FRET_distance_restraint'] = \
                cur_fret_distance_restraint_index

        # => FPS AV parameters
        # donor
        for i in range(nr_entries):
            # only do this if there are the linker length, the linker width and at least one probe radius given.
            if (('FLR_FPS_AV_param_donor_AV_linker_length' in xls_data.keys()
                 and not pandas.isnull(xls_data['FLR_FPS_AV_param_donor_AV_linker_length'][i]))
                    and ('FLR_FPS_AV_param_donor_AV_linker_width' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_AV_param_donor_AV_linker_width'][i]))
                    and ('FLR_FPS_AV_param_donor_AV_probe_radius_1' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_AV_param_donor_AV_probe_radius_1'][i]))):
                cur_FPS_AV_param_donor_AV_num_linker_atoms = None if (
                        'FLR_FPS_AV_param_donor_AV_num_linker_atoms' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_FPS_AV_param_donor_AV_num_linker_atoms'][i])) \
                    else int(xls_data['FLR_FPS_AV_param_donor_AV_num_linker_atoms'][i])
                cur_FPS_AV_param_donor_AV_linker_length = \
                    float(xls_data['FLR_FPS_AV_param_donor_AV_linker_length'][i])
                cur_FPS_AV_param_donor_AV_linker_width = \
                    float(xls_data['FLR_FPS_AV_param_donor_AV_linker_width'][i])
                cur_FPS_AV_param_donor_AV_probe_radius_1 = \
                    float(xls_data['FLR_FPS_AV_param_donor_AV_probe_radius_1'][i])
                cur_FPS_AV_param_donor_AV_probe_radius_2 = None if (
                        'FLR_FPS_AV_param_donor_AV_probe_radius_2' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_FPS_AV_param_donor_AV_probe_radius_2'][i])) \
                    else float(xls_data['FLR_FPS_AV_param_donor_AV_probe_radius_2'][i])
                cur_FPS_AV_param_donor_AV_probe_radius_3 = None if (
                        'FLR_FPS_AV_param_donor_AV_probe_radius_3' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_FPS_AV_param_donor_AV_probe_radius_3'][i])) \
                    else float(xls_data['FLR_FPS_AV_param_donor_AV_probe_radius_3'][i])
                cur_FPS_AV_param_donor = ihm.flr.FPSAVParameter(
                    num_linker_atoms=cur_FPS_AV_param_donor_AV_num_linker_atoms,
                    linker_length=cur_FPS_AV_param_donor_AV_linker_length,
                    linker_width=cur_FPS_AV_param_donor_AV_linker_width,
                    probe_radius_1=cur_FPS_AV_param_donor_AV_probe_radius_1,
                    probe_radius_2=cur_FPS_AV_param_donor_AV_probe_radius_2,
                    probe_radius_3=cur_FPS_AV_param_donor_AV_probe_radius_3)

                cur_FPS_AV_param_donor_index = -1
                if cur_FPS_AV_param_donor not in self.list_FPS_AV_parameters:
                    self.list_FPS_AV_parameters.append(cur_FPS_AV_param_donor)
                cur_FPS_AV_param_donor_index = \
                    self.list_FPS_AV_parameters.index(cur_FPS_AV_param_donor)
                self.list_of_object_indices[i]['FPS_AV_param_donor'] = \
                    cur_FPS_AV_param_donor_index
                # save whether AV1 or AV3
                self.list_of_object_indices[i]['AV_modeling_method_donor'] = 'AV1' \
                    if (cur_FPS_AV_param_donor_AV_probe_radius_2 == None) else 'AV3'
        # acceptor
        for i in range(nr_entries):
            # only do this if there are the linker length, the linker width
            # and at least one probe radius given.
            if (('FLR_FPS_AV_param_acceptor_AV_linker_length' in xls_data.keys()
                 and not pandas.isnull(xls_data['FLR_FPS_AV_param_acceptor_AV_linker_length'][i]))
                    and ('FLR_FPS_AV_param_acceptor_AV_linker_width' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_AV_param_acceptor_AV_linker_width'][i]))
                    and ('FLR_FPS_AV_param_acceptor_AV_probe_radius_1' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_AV_param_acceptor_AV_probe_radius_1'][i]))):
                cur_FPS_AV_param_acceptor_AV_num_linker_atoms = None if (
                        'FLR_FPS_AV_param_acceptor_AV_num_linker_atoms' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_FPS_AV_param_acceptor_AV_num_linker_atoms'][i])) \
                    else int(xls_data['FLR_FPS_AV_param_acceptor_AV_num_linker_atoms'][i])
                cur_FPS_AV_param_acceptor_AV_linker_length = \
                    float(xls_data['FLR_FPS_AV_param_acceptor_AV_linker_length'][i])
                cur_FPS_AV_param_acceptor_AV_linker_width = \
                    float(xls_data['FLR_FPS_AV_param_acceptor_AV_linker_width'][i])
                cur_FPS_AV_param_acceptor_AV_probe_radius_1 = \
                    float(xls_data['FLR_FPS_AV_param_acceptor_AV_probe_radius_1'][i])
                cur_FPS_AV_param_acceptor_AV_probe_radius_2 = None if (
                        'FLR_FPS_AV_param_acceptor_AV_probe_radius_2' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_FPS_AV_param_acceptor_AV_probe_radius_2'][i])) \
                    else float(xls_data['FLR_FPS_AV_param_acceptor_AV_probe_radius_2'][i])
                cur_FPS_AV_param_acceptor_AV_probe_radius_3 = None if (
                        'FLR_FPS_AV_param_acceptor_AV_probe_radius_3' not in xls_data.keys()
                        or pandas.isnull(xls_data['FLR_FPS_AV_param_acceptor_AV_probe_radius_3'][i])) \
                    else float(xls_data['FLR_FPS_AV_param_acceptor_AV_probe_radius_3'][i])
                cur_FPS_AV_param_acceptor = ihm.flr.FPSAVParameter(
                    num_linker_atoms=cur_FPS_AV_param_acceptor_AV_num_linker_atoms,
                    linker_length=cur_FPS_AV_param_acceptor_AV_linker_length,
                    linker_width=cur_FPS_AV_param_acceptor_AV_linker_width,
                    probe_radius_1=cur_FPS_AV_param_acceptor_AV_probe_radius_1,
                    probe_radius_2=cur_FPS_AV_param_acceptor_AV_probe_radius_2,
                    probe_radius_3=cur_FPS_AV_param_acceptor_AV_probe_radius_3)

                cur_FPS_AV_param_acceptor_index = -1
                if cur_FPS_AV_param_acceptor not in self.list_FPS_AV_parameters:
                    self.list_FPS_AV_parameters.append(cur_FPS_AV_param_acceptor)
                cur_FPS_AV_param_acceptor_index = \
                    self.list_FPS_AV_parameters.index(cur_FPS_AV_param_acceptor)
                self.list_of_object_indices[i]['FPS_AV_param_acceptor'] = \
                    cur_FPS_AV_param_acceptor_index
                ## save whether it is AV1 or AV3
                self.list_of_object_indices[i]['AV_modeling_method_acceptor'] = 'AV1' if (
                        cur_FPS_AV_param_acceptor_AV_probe_radius_2 == None) else 'AV3'

        # => AV Modeling

        # ==> FPS_modeling_AV
        for i in range(nr_entries):
            # for each modeling step that include FPS
            for this_ihm_modeling_protocol_id in self.list_ihm_modeling_protocols:
                for this_step in this_ihm_modeling_protocol_id.steps:
                    if 'FPS' in this_step.software.name or 'Olga' in this_step.software.name:
                        # The fps global parameters are stored by protocol id
                        # and then by step id. To get to this, we have to go
                        # via the protocol id that is stored in
                        # self.list_ihm_modeling_protocols_ids
                        this_fps_global_parameters = \
                            self.dict_flr_fps_global_parameters_by_protocol_id[
                                self.list_ihm_modeling_protocols_ids[
                                    self.list_ihm_modeling_protocols.index(
                                        this_ihm_modeling_protocol_id)]][
                                self.list_ihm_modeling_protocol_modeling_step_ids[
                                    self.list_ihm_modeling_protocol_modeling_steps.index(
                                        this_step)]]
                        # donor
                        if 'AV_modeling_method_donor' in self.list_of_object_indices[i].keys():
                            cur_FPS_modeling_AV = ihm.flr.FPSModeling(
                                protocol=this_step,
                                restraint_group=self.list_fret_distance_restraint_groups[
                                    self.list_of_object_indices[i]['FRET_distance_restraint_group']],
                                global_parameter=this_fps_global_parameters,
                                probe_modeling_method=self.list_of_object_indices[i]['AV_modeling_method_donor'],
                                details=this_step.name)
                            cur_FPS_modeling_AV_index = -1
                            if cur_FPS_modeling_AV not in self.list_FPS_modeling_by_AV:
                                self.list_FPS_modeling_by_AV.append(cur_FPS_modeling_AV)
                            cur_FPS_modeling_AV_index = \
                                self.list_FPS_modeling_by_AV.index(cur_FPS_modeling_AV)
                            if 'FPS_modeling_AV_donor' not in self.list_of_object_indices[i].keys():
                                self.list_of_object_indices[i]['FPS_modeling_AV_donor'] = []
                            self.list_of_object_indices[i][
                                'FPS_modeling_AV_donor'].append(cur_FPS_modeling_AV_index)

                        # acceptor
                        if 'AV_modeling_method_acceptor' in self.list_of_object_indices[i].keys():
                            cur_FPS_modeling_AV = ihm.flr.FPSModeling(
                                protocol=this_step,
                                restraint_group=self.list_fret_distance_restraint_groups[
                                    self.list_of_object_indices[i]['FRET_distance_restraint_group']],
                                global_parameter=this_fps_global_parameters,
                                probe_modeling_method=self.list_of_object_indices[i]['AV_modeling_method_acceptor'],
                                details=this_step.name)
                            cur_FPS_modeling_AV_index = -1
                            if cur_FPS_modeling_AV not in self.list_FPS_modeling_by_AV:
                                self.list_FPS_modeling_by_AV.append(cur_FPS_modeling_AV)
                            cur_FPS_modeling_AV_index = \
                                self.list_FPS_modeling_by_AV.index(cur_FPS_modeling_AV)
                            if 'FPS_modeling_AV_acceptor' not in self.list_of_object_indices[i].keys():
                                self.list_of_object_indices[i]['FPS_modeling_AV_acceptor'] = []
                            self.list_of_object_indices[i][
                                'FPS_modeling_AV_acceptor'].append(cur_FPS_modeling_AV_index)

        for i in range(nr_entries):
            # donor
            if 'FPS_modeling_AV_donor' in self.list_of_object_indices[i].keys():
                for this_FPS_modeling_AV_donor_entry in \
                        self.list_of_object_indices[i]['FPS_modeling_AV_donor']:
                    cur_FPS_AV_modeling_FPS_modeling_id = \
                        self.list_FPS_modeling_by_AV[this_FPS_modeling_AV_donor_entry]
                    cur_FPS_AV_modeling_sample_probe_id_donor = \
                        self.list_sample_probe_details[
                            self.list_of_object_indices[i]['Sample_probe_details_donor']]
                    cur_FPS_AV_modeling_parameter_id_donor = \
                        self.list_FPS_AV_parameters[
                            self.list_of_object_indices[i]['FPS_AV_param_donor']]

                    cur_FPS_AV_modeling_donor = ihm.flr.FPSAVModeling(
                        fps_modeling=cur_FPS_AV_modeling_FPS_modeling_id,
                        sample_probe=cur_FPS_AV_modeling_sample_probe_id_donor,
                        parameter=cur_FPS_AV_modeling_parameter_id_donor)
                    cur_FPS_AV_modeling_donor_index = -1
                    if cur_FPS_AV_modeling_donor not in self.list_FPS_AV_modeling:
                        self.list_FPS_AV_modeling.append(cur_FPS_AV_modeling_donor)
                    cur_FPS_AV_modeling_donor_index = \
                        self.list_FPS_AV_modeling.index(cur_FPS_AV_modeling_donor)
                    self.list_of_object_indices[i]['FPS_AV_modeling_AV_donor'] = \
                        cur_FPS_AV_modeling_donor_index

            # acceptor
            if 'FPS_modeling_AV_acceptor' in self.list_of_object_indices[i].keys():
                for this_FPS_modeling_AV_acceptor_entry in \
                        self.list_of_object_indices[i]['FPS_modeling_AV_acceptor']:
                    cur_FPS_AV_modeling_FPS_modeling_id = \
                        self.list_FPS_modeling_by_AV[this_FPS_modeling_AV_acceptor_entry]
                    cur_FPS_AV_modeling_sample_probe_id_acceptor = \
                        self.list_sample_probe_details[
                            self.list_of_object_indices[i]['Sample_probe_details_acceptor']]
                    cur_FPS_AV_modeling_parameter_id_acceptor = \
                        self.list_FPS_AV_parameters[
                            self.list_of_object_indices[i]['FPS_AV_param_acceptor']]
                    cur_FPS_AV_modeling_acceptor = ihm.flr.FPSAVModeling(
                        fps_modeling=cur_FPS_AV_modeling_FPS_modeling_id,
                        sample_probe=cur_FPS_AV_modeling_sample_probe_id_acceptor,
                        parameter=cur_FPS_AV_modeling_parameter_id_acceptor)
                    cur_FPS_AV_modeling_acceptor_index = -1
                    if cur_FPS_AV_modeling_acceptor not in self.list_FPS_AV_modeling:
                        self.list_FPS_AV_modeling.append(cur_FPS_AV_modeling_acceptor)
                    cur_FPS_AV_modeling_acceptor_index = \
                        self.list_FPS_AV_modeling.index(cur_FPS_AV_modeling_acceptor)
                    self.list_of_object_indices[i]['FPS_AV_modeling_AV_acceptor'] = \
                        cur_FPS_AV_modeling_acceptor_index

        # ==> FPS_mean_probe_position
        # Get the mpp reference group and the coordinates of the position
        for i in range(nr_entries):
            # donor
            # Check whether all columns are present and not empty
            if (('FLR_FPS_MPP_reference_group_id_donor' in xls_data.keys()
                 and not pandas.isnull(xls_data['FLR_FPS_MPP_reference_group_id_donor'][i]))
                    and ('FLR_FPS_MPP_xcoord_donor' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_MPP_xcoord_donor'][i]))
                    and ('FLR_FPS_MPP_ycoord_donor' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_MPP_ycoord_donor'][i]))
                    and ('FLR_FPS_MPP_zcoord_donor' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_MPP_zcoord_donor'][i]))):
                cur_FPS_MPP_reference_group_id = xls_data['FLR_FPS_MPP_reference_group_id_donor'][i]
                cur_FPS_MPP_xcoord = xls_data['FLR_FPS_MPP_xcoord_donor'][i]
                cur_FPS_MPP_ycoord = xls_data['FLR_FPS_MPP_ycoord_donor'][i]
                cur_FPS_MPP_zcoord = xls_data['FLR_FPS_MPP_zcoord_donor'][i]
                cur_FPS_MPP_donor = ihm.flr.FPSMeanProbePosition(
                    sample_probe=self.list_sample_probe_details[
                        self.list_of_object_indices[i]['Sample_probe_details_donor']],
                    x=cur_FPS_MPP_xcoord,
                    y=cur_FPS_MPP_ycoord,
                    z=cur_FPS_MPP_zcoord)

                cur_FPS_mean_probe_position_donor_index = -1
                if cur_FPS_MPP_donor not in self.list_FPS_mean_probe_positions:
                    self.list_FPS_mean_probe_positions.append(cur_FPS_MPP_donor)
                cur_FPS_mean_probe_position_donor_index = \
                    self.list_FPS_mean_probe_positions.index(cur_FPS_MPP_donor)
                self.list_of_object_indices[i]['FPS_mean_probe_position_donor'] = \
                    cur_FPS_mean_probe_position_donor_index

            # acceptor
            # Check whether all columns are present and not empty
            if (('FLR_FPS_MPP_reference_group_id_acceptor' in xls_data.keys()
                 and not pandas.isnull(xls_data['FLR_FPS_MPP_reference_group_id_acceptor'][i]))
                    and ('FLR_FPS_MPP_xcoord_acceptor' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_MPP_xcoord_acceptor'][i]))
                    and ('FLR_FPS_MPP_ycoord_acceptor' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_MPP_ycoord_acceptor'][i]))
                    and ('FLR_FPS_MPP_zcoord_acceptor' in xls_data.keys()
                         and not pandas.isnull(xls_data['FLR_FPS_MPP_zcoord_acceptor'][i]))):
                cur_FPS_MPP_reference_group_id = xls_data['FLR_FPS_MPP_reference_group_id_acceptor'][i]
                cur_FPS_MPP_xcoord = xls_data['FLR_FPS_MPP_xcoord_acceptor'][i]
                cur_FPS_MPP_ycoord = xls_data['FLR_FPS_MPP_ycoord_acceptor'][i]
                cur_FPS_MPP_zcoord = xls_data['FLR_FPS_MPP_zcoord_acceptor'][i]
                cur_FPS_MPP_acceptor = ihm.flr.FPSMeanProbePosition(
                    sample_probe=self.list_sample_probe_details[
                        self.list_of_object_indices[i]['Sample_probe_details_acceptor']],
                    x=cur_FPS_MPP_xcoord,
                    y=cur_FPS_MPP_ycoord,
                    z=cur_FPS_MPP_zcoord)

                cur_FPS_mean_probe_position_acceptor_index = -1
                if cur_FPS_MPP_acceptor not in self.list_FPS_mean_probe_positions:
                    self.list_FPS_mean_probe_positions.append(cur_FPS_MPP_acceptor)
                cur_FPS_mean_probe_position_acceptor_index = \
                    self.list_FPS_mean_probe_positions.index(cur_FPS_MPP_acceptor)
                self.list_of_object_indices[i]['FPS_mean_probe_position_acceptor'] = \
                    cur_FPS_mean_probe_position_acceptor_index

        # ==> FPS_modeling by MPP
        for i in range(nr_entries):
            # For each modeling step that include FPS
            for this_ihm_modeling_protocol_id in self.list_ihm_modeling_protocols:
                for this_step in this_ihm_modeling_protocol_id.steps:
                    if 'FPS' in this_step.software.name or 'Olga' in this_step.software.name:
                        # the fps global parameters are stored by protocol id
                        # and then by step id
                        # to get to this, we have to go via the protocol id
                        # that is stored in list_ihm_modeling_protocols_ids
                        this_fps_global_parameters = \
                            self.dict_flr_fps_global_parameters_by_protocol_id[
                                self.list_ihm_modeling_protocols_ids[
                                    self.list_ihm_modeling_protocols.index(
                                        this_ihm_modeling_protocol_id)]][
                                self.list_ihm_modeling_protocol_modeling_step_ids[
                                    self.list_ihm_modeling_protocol_modeling_steps.index(
                                        this_step)]]
                        # donor
                        if 'FPS_mean_probe_position_donor' in self.list_of_object_indices[i].keys():
                            cur_FPS_modeling_MPP = ihm.flr.FPSModeling(
                                protocol=this_step,
                                restraint_group=
                                self.list_fret_distance_restraint_groups[
                                    self.list_of_object_indices[i]['FRET_distance_restraint_group']],
                                global_parameter=this_fps_global_parameters,
                                probe_modeling_method='MPP',
                                details=this_step.name)
                            cur_FPS_modeling_MPP_index = -1
                            if cur_FPS_modeling_MPP not in self.list_FPS_modeling_by_MPP:
                                self.list_FPS_modeling_by_MPP.append(cur_FPS_modeling_MPP)
                            cur_FPS_modeling_MPP_index = \
                                self.list_FPS_modeling_by_MPP.index(cur_FPS_modeling_MPP)
                            if 'FPS_modeling_MPP_donor' not in self.list_of_object_indices[i].keys():
                                self.list_of_object_indices[i]['FPS_modeling_MPP_donor'] = []
                            self.list_of_object_indices[i][
                                'FPS_modeling_MPP_donor'].append(cur_FPS_modeling_MPP_index)
                        # acceptor
                        if 'FPS_mean_probe_position_acceptor' in self.list_of_object_indices[i].keys():
                            cur_FPS_modeling_MPP = ihm.flr.FPSModeling(
                                protocol=this_step,
                                restraint_group=
                                self.list_fret_distance_restraint_groups[
                                    self.list_of_object_indices[i]['FRET_distance_restraint_group']],
                                global_parameter=this_fps_global_parameters,
                                probe_modeling_method='MPP',
                                details=this_step.name)
                            cur_FPS_modeling_MPP_index = -1
                            if cur_FPS_modeling_MPP not in self.list_FPS_modeling_by_MPP:
                                self.list_FPS_modeling_by_MPP.append(cur_FPS_modeling_MPP)
                            cur_FPS_modeling_MPP_index = \
                                self.list_FPS_modeling_by_MPP.index(cur_FPS_modeling_MPP)
                            if 'FPS_modeling_MPP_acceptor' not in self.list_of_object_indices[i].keys():
                                self.list_of_object_indices[i]['FPS_modeling_MPP_acceptor'] = []
                            self.list_of_object_indices[i][
                                'FPS_modeling_MPP_acceptor'].append(cur_FPS_modeling_MPP_index)
        # Connect FPS_modeling with the Mean probe position and the mpp_atom_position_group
        for i in range(nr_entries):
            # donor
            if 'FPS_modeling_MPP_donor' in self.list_of_object_indices[i].keys():
                for this_FPS_modeling_MPP_donor_entry in \
                        self.list_of_object_indices[i]['FPS_modeling_MPP_donor']:
                    cur_FPS_MPP_reference_group_id = \
                        xls_data['FLR_FPS_MPP_reference_group_id_donor'][i]
                    cur_FPS_MPP_modeling_FPS_modeling_id = \
                        self.list_FPS_modeling_by_MPP[this_FPS_modeling_MPP_donor_entry]
                    cur_FPS_MPP_modeling_MPP_id = \
                        self.list_FPS_mean_probe_positions[
                            self.list_of_object_indices[i]['FPS_mean_probe_position_donor']]
                    cur_FPS_MPP_modeling_MPP_group_id = \
                        self.list_flr_fps_mpp_groups[cur_FPS_MPP_reference_group_id]
                    cur_FPS_MPP_modeling = ihm.flr.FPSMPPModeling(
                        fps_modeling=cur_FPS_MPP_modeling_FPS_modeling_id,
                        mpp=cur_FPS_MPP_modeling_MPP_id,
                        mpp_atom_position_group=cur_FPS_MPP_modeling_MPP_group_id)
                    cur_FPS_MPP_modeling_donor_index = -1
                    if cur_FPS_MPP_modeling not in self.list_FPS_MPP_modeling:
                        self.list_FPS_MPP_modeling.append(cur_FPS_MPP_modeling)
                    cur_FPS_MPP_modeling_donor_index = \
                        self.list_FPS_MPP_modeling.index(cur_FPS_MPP_modeling)
                    self.list_of_object_indices[i]['FPS_MPP_modeling_donor'] = \
                        cur_FPS_MPP_modeling_donor_index

            # acceptor
            if 'FPS_modeling_MPP_acceptor' in self.list_of_object_indices[i].keys():
                for this_FPS_modeling_MPP_acceptor_entry in \
                        self.list_of_object_indices[i]['FPS_modeling_MPP_acceptor']:
                    cur_FPS_MPP_reference_group_id = xls_data['FLR_FPS_MPP_reference_group_id_acceptor'][i]
                    cur_FPS_MPP_modeling_FPS_modeling_id = \
                        self.list_FPS_modeling_by_MPP[this_FPS_modeling_MPP_acceptor_entry]
                    cur_FPS_MPP_modeling_MPP_id = \
                        self.list_FPS_mean_probe_positions[
                            self.list_of_object_indices[i]['FPS_mean_probe_position_acceptor']]
                    cur_FPS_MPP_modeling_MPP_group_id = \
                        self.list_flr_fps_mpp_groups[cur_FPS_MPP_reference_group_id]
                    cur_FPS_MPP_modeling = ihm.flr.FPSMPPModeling(
                        fps_modeling=cur_FPS_MPP_modeling_FPS_modeling_id,
                        mpp=cur_FPS_MPP_modeling_MPP_id,
                        mpp_atom_position_group=cur_FPS_MPP_modeling_MPP_group_id)
                    cur_FPS_MPP_modeling_acceptor_index = -1
                    if cur_FPS_MPP_modeling not in self.list_FPS_MPP_modeling:
                        self.list_FPS_MPP_modeling.append(cur_FPS_MPP_modeling)
                    cur_FPS_MPP_modeling_acceptor_index = \
                        self.list_FPS_MPP_modeling.index(cur_FPS_MPP_modeling)
                    self.list_of_object_indices[i]['FPS_MPP_modeling_acceptor'] = \
                        cur_FPS_MPP_modeling_acceptor_index


    def add_flr_model_quality(self,
                              sheet_name='FLR_FRET_Model_quality',
                              skiprows=3,
                              header=0):
        """Read information on the FLR model quality from the excel sheet"""
        if self.verbose:
            print(" ... Processing tab \'FLR_FRET_Model_quality\' ...")
        xls_data = pandas.pandas.read_excel(
            self.xls_file,
            sheet_name=sheet_name,
            skiprows=skiprows,
            header=header)
        nr_entries = len(xls_data['FLR_FRET_model_quality_Ordinal_id'])

        for i in range(nr_entries):
            cur_flr_model_quality_ordinal_id = \
                xls_data['FLR_FRET_model_quality_Ordinal_id'][i]
            cur_flr_model_quality_model_id = xls_data['FLR_FRET_model_quality_Model_id'][i]
            cur_flr_model_quality_chi_square_reduced = \
                xls_data['FLR_FRET_model_quality_Chi_square_reduced'][i]
            cur_flr_model_quality_dataset_group_id = \
            xls_data['FLR_FRET_model_quality_dataset_group_id'][i]
            cur_flr_model_quality_method = None if (
                    'FLR_FRET_model_quality_Method' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_model_quality_Method'][i])) \
                else xls_data['FLR_FRET_model_quality_Method'][i]
            cur_flr_model_quality_details = None if (
                    'FLR_FRET_model_quality_details' not in xls_data.keys()
                    or pandas.isnull(xls_data['FLR_FRET_model_quality_details'][i])) \
                else xls_data['FLR_FRET_model_quality_details'][i]

            cur_flr_model_quality = ihm.flr.FRETModelQuality(
                model=self.list_models[
                    self.list_models_ids.index(cur_flr_model_quality_model_id)],
                chi_square_reduced=cur_flr_model_quality_chi_square_reduced,
                dataset_group=self.list_dataset_groups[
                    self.list_dataset_group_ids.index(
                        cur_flr_model_quality_dataset_group_id)],
                method=cur_flr_model_quality_method,
                details=cur_flr_model_quality_details)
            if cur_flr_model_quality not in self.list_flr_model_quality:
                self.list_flr_model_quality.append(cur_flr_model_quality)

    def add_flr_model_distances(self,
                                sheet_name='FLR_FRET_Model_distances',
                                skiprows=3,
                                header=0):
        """Read the FLR model distances from the excel sheet."""
        if self.verbose:
            print(' ... Processing tab \'FLR_FRET_Model_distances\' ...')
        xls_data = pandas.read_excel(self.xls_file,
                                     sheet_name=sheet_name,
                                     skiprows=skiprows,
                                     header=header)
        nr_entries = len(xls_data['FLR_FRET_model_distance_Ordinal_id'])
        for i in range(nr_entries):
            cur_flr_model_distance_ordinal_id = \
                xls_data['FLR_FRET_model_distance_Ordinal_id'][i]
            cur_flr_model_distance_model_id = \
                xls_data['FLR_FRET_model_distance_Model_id'][i]
            cur_flr_model_distance_restraint_id = \
                xls_data['FLR_FRET_model_distance_Distance_restraint_id'][i]
            cur_flr_model_distance_distance = \
                xls_data['FLR_FRET_model_distance_Model_distance'][i]

            cur_flr_fret_model_distance = ihm.flr.FRETModelDistance(
                restraint=self.list_fret_distance_restraints[
                    self.list_fret_distance_restraint_ids.index(
                        cur_flr_model_distance_restraint_id)],
                model=self.list_models[self.list_models_ids.index(
                    cur_flr_model_distance_model_id)],
                distance=cur_flr_model_distance_distance)

            if not cur_flr_fret_model_distance in self.list_flr_model_distances:
                self.list_flr_model_distances.append(cur_flr_fret_model_distance)

    def create_flrdata(self):
        """Create the FLRData object"""
        FLR_collection1 = ihm.flr.FLRData()
        # add FRET distance restraints
        for entry in self.list_fret_distance_restraint_groups:
            FLR_collection1.distance_restraint_groups.append(entry)

        # add the poly_probe_conjugates
        for entry in self.list_poly_probe_conjugates:
            FLR_collection1.poly_probe_conjugates.append(entry)

        # add FPS modeling
        for entry in self.list_FPS_AV_modeling:
            FLR_collection1.fps_modeling.append(entry)
        for entry in self.list_FPS_MPP_modeling:
            FLR_collection1.fps_modeling.append(entry)

        # add the FRET model quality
        for entry in self.list_flr_model_quality:
            FLR_collection1.fret_model_qualities.append(entry)

        # add the FRET model distances
        for entry in self.list_flr_model_distances:
            FLR_collection1.fret_model_distances.append(entry)

        self.system.flr_data = [FLR_collection1]

        # Check for struct_assemblies that were not used elsewhere and add
        # them to the system.orphan_assemblies
        seen_assemblies = []
        for entry in self.list_ihm_modeling_protocol_analysis_steps:
            seen_assemblies.append(entry.assembly)
        for entry in self.list_ihm_modeling_protocol_modeling_steps:
            seen_assemblies.append(entry.assembly)
        for entry in self.list_structure_assemblies:
            if not occurs_in_list(entry, seen_assemblies):
                self.system.orphan_assemblies.append(entry)

    def write_cif_file(self):
        """Write the system to a cif file."""
        if self.verbose:
            print('>>> Writing \'%s\'' % (self.cifout_filename))
        with open(self.cifout_filename, 'w') as fh:
            ihm.dumper.write(fh, [self.system])

        if self.verbose:
            print('Done.')


    def handle_atom_site_file(self):
        """Handle the atom_site entries for the models.
        It is important to add the atoms to the models, since models in the
        output might not be in the order as they are read.
        This is due to the models being collected e.g. by state.
        """
        # First read the mmcif file
        # This file has to contain _ihm_model_list.model_id,
        # _ihm_model_group.id, and _ihm_model_group_link.group_id
        # together with _ihm_model_group_link.model_id.
        # The _ihm_model_group.id can be only one, but all model_ids have to
        # be assigned to this model_group.
        # This collection of models is separate from the system.
        # Only the atoms for the models will be taken from it.

        # First check whether all the information is in the atom_site file.
        # If not everything requires is present, it is added to a new file
        new_atom_site_filename = self.check_atom_site_file_for_model_list(
            self.atom_site_filename)
        self.atom_site_filename = new_atom_site_filename
        if self.verbose:
            print('<<< Adding _atom_site information from file \'%s\''
                  % (self.atom_site_filename))
            print('Note: It is assumed that the order of the models '
                  'corresponds to the model IDs given in the excel sheet!')
        temp_atom_site_systems = ihm.System()
        with open(self.atom_site_filename) as fh:
            temp_atom_site_systems = ihm.reader.read(fh)
        cur_temp_atom_site_system = temp_atom_site_systems[0]
        # We assume that the model_ids are in the same order as the model
        # ids that are given int he excel sheet.
        # system._all_models() returns model_group and model.
        # So we only take the second entry in the tuple
        list_models_ids_string = [str(x) for x in self.list_models_ids]
        for this_model_group, this_model in \
                cur_temp_atom_site_system._all_models():
            # Get the model_id that is assigned to the model when it is read.
            this_model_id = this_model._id
            # now find the respective model
            # !!! Note: Converting the id here to an int assumes that the
            # model_ids in the Model tab of the excel sheet are int as well.
            # Therefore, convert all the ids to strings before and now
            # compare strings with strings.
            index_of_this_model_id = \
                list_models_ids_string.index(str(this_model_id))

            cur_model = self.list_models[index_of_this_model_id]
            # Now get the atoms from this_model and add them to cur_model
            for this_atom in this_model.get_atoms():
                for this_asym_unit in self.list_asym_units:
                    # sometimes this_atom.asym_unit.id is set and sometimes
                    # it is not. If it is not set, then
                    # this_atom.asym_unit._id is set
                    if this_atom.asym_unit.id is None:
                        if this_atom.asym_unit._id == this_asym_unit.id:
                            this_atom.asym_unit = this_asym_unit
                            break
                    else:
                        if this_atom.asym_unit.id == this_asym_unit.id:
                            this_atom.asym_unit = this_asym_unit
                            break
                cur_model.add_atom(this_atom)

    def read_atom_site_entries(self,atom_site_file):
        """ Read the atom site entries from a mmcif file. """
        # read all the data
        print("<<< Reading \"%s\"" % (atom_site_file))
        data = []
        infile = open(atom_site_file, 'r')
        for line in infile.readlines():
            data.append(line)
        infile.close()
        # sort out the entries that do not belong to the atom_site section
        atom_site_data = []
        atom_site_found = False
        for entry in data:
            # if the atom_site entries were encountered
            if atom_site_found:
                # if encounter the end of the atom_site entries
                if '#' in entry:
                    break
                else:
                    atom_site_data.append(entry)
                continue
            else:
                if '_atom_site.' in entry:
                    atom_site_found = True
                    atom_site_data.append(entry)
                    continue
                else:
                    continue
        return atom_site_data

    def do_add_atom_site_entries(self, atom_site_file, mmcif_file):
        """ Reads the atom site entries from a mmcif file and adds it to the
          previously created mmcif file. Additionally, it adds the entry
          "ihm_model_id" with the value of pdbx_PDB_model_num.
          ! Note: The order of the models here has to be the same as used
          in the Excel sheet !
        """
        print(
            '... Adding atom_site entries ... Note: This step will only '
            'extract the atom_site entries and omit all other information.')
        # Read the entries from the atom_site_file
        read_data = self.read_atom_site_entries(atom_site_file)
        # split header and ATOM entries
        header = []
        atom_entries = []
        for entry in read_data:
            # From the reading routine, there should either be the header of
            # the section or the atom entries
            if '_atom_site' in entry:
                header.append(entry)
            else:
                atom_entries.append(entry)
        # add a row to the header _atom_site.ihm_model_id
        header.append('_atom_site.ihm_model_id\n')
        # modify each atom entry by adding the ihm_model_id that corresponds
        # to the pdbx_PDB_model_num
        # get the index of the pdbx_PDB_model_num
        index_pdbx_PDB_model_num = \
            [header.index(entry) for entry in header
             if '_atom_site.pdbx_PDB_model_num' in entry]
        if index_pdbx_PDB_model_num == []:
            print(
                'ERROR: The file \'%s\' does not contain '
                '_atom_site.pdbx_PDB_model_num. Please add this information '
                'to your file. Atom_site entries will not be written.\n')
            return
        index_pdbx_PDB_model_num = index_pdbx_PDB_model_num[0]
        # get the values from as _atom_site.pdbx_PDB_model_num
        # and add a column to the data that has the same value as
        # _atom_site.pdbx_PDB_model_num
        new_atom_entries = []
        for entry in atom_entries:
            splitentry = entry.split()
            cur_pdbx_pdb_model_num = splitentry[index_pdbx_PDB_model_num]
            # remove the endline character and add the new entry
            new_entry = '%s\t%s\n' % (entry.split('\n')[0],
                                      cur_pdbx_pdb_model_num)
            new_atom_entries.append(new_entry)
        # Collect atom types
        index_type_symbol = [header.index(entry) for entry in header if
                             '_atom_site.type_symbol' in entry]
        if index_type_symbol == []:
            print(
                'ERROR: The file \'%s\' does not contain '
                '_atom_site.type_symbol. Please add this information to your '
                'file. Atom_site entries will not be written.\n')
            return
        index_type_symbol = index_type_symbol[0]
        atom_types = []
        for entry in atom_entries:
            splitentry = entry.split()
            cur_type_symbol = splitentry[index_type_symbol]
            if cur_type_symbol not in atom_types:
                atom_types.append(cur_type_symbol)

        # write the header and the new entries to the mmcif_file
        print('>>> Writing atom_site entries to \'%s\'' % (mmcif_file))
        outfile = open(mmcif_file, 'a')
        # Write the atom types
        outfile.write('#\n')
        outfile.write('loop_\n')
        outfile.write('_atom_type.symbol\n')
        for entry in atom_types:
            outfile.write('%s\n' % (entry))
        outfile.write('#\n')
        ## Write the atom_site entries
        outfile.write('#\n')
        outfile.write('loop_\n')
        for entry in header:
            outfile.write(entry)
        for entry in new_atom_entries:
            outfile.write(entry)
        outfile.write('#\n')
        outfile.close()

    def check_atom_site_file_for_model_list(self, atom_site_filename):
        """
        Check the file with the atom_site entries, whether it includes
        information on _ihm_model_list, _ihm_model_group, and
        _ihm_model_group_link entries.
        This is important, because python-ihm will not read the atom_site
        entries otherwise. If it is not present, then add the information
        and write a new file.
        :param atom_site_filename: Filename of the mmcif file containing the
        atom_site entries.
        :return: Filename of the mmcif file to use for the atom_site entries.
        """
        if self.verbose:
            print(".. Checking atom_site file \'%s\' .." % (
                atom_site_filename))
        infile = open(atom_site_filename, 'r')
        # Were all the three necessary categories found?
        ihm_model_list_found = False
        ihm_model_group_found = False
        ihm_model_group_link_found = False
        # Collect the header entries for the atom_site entries
        atom_site_header_entries = []
        atom_site_entries_found = False
        atom_site_entries = []
        for line in infile.readlines():
            if "_ihm_model_list" in line:
                ihm_model_list_found = True
            if "_ihm_model_group" in line:
                ihm_model_group_found = True
            if "_ihm_model_group_link" in line:
                ihm_model_group_link_found = True
            if "_atom_site" in line:
                atom_site_header_entries.append(line.split()[0])
                atom_site_entries_found = True
            # if _atom_site header was already encountered, store the
            # atom_site entries
            if atom_site_entries_found:
                # The end of the atom_site_entries
                if "#" in line:
                    atom_site_entries_found = False
                if "ATOM" in line:
                    atom_site_entries.append(line)

        # Alternatively here: omit the given information on the models
        # If all three categories were found
        if ihm_model_list_found \
                and ihm_model_group_found \
                and ihm_model_group_link_found:
            # TODO: Check whether the entries are correct and complete
            # For now, assume that all the entries are complete. So, simply
            # return the original filename
            return atom_site_filename
        # Otherwise, keep the atom_site entries and count the number of models
        else:
            # get the pdbx_PDB_model_num index to find the number of models
            pdbx_pdb_model_num_index = \
                atom_site_header_entries.index("_atom_site.pdbx_PDB_model_num")
            # now go through the atom_site entries to collect the model numbers
            model_numbers = []
            for entry in atom_site_entries:
                # get the current model number
                cur_model_number = entry.split()[pdbx_pdb_model_num_index]
                if cur_model_number not in model_numbers:
                    model_numbers.append(cur_model_number)
            number_of_models = len(model_numbers)
            new_lines = []
            # _ihm_model_list
            new_lines.append("loop_\n")
            new_lines.append("_ihm_model_list.model_id\n")
            for entry in model_numbers:
                new_lines.append("%s\n" % (entry))
            new_lines.append("#\n#\n")
            # _ihm_model_group (only one model_group)
            new_lines.append("loop_\n")
            new_lines.append("_ihm_model_group.id\n1\n#\n#\n")
            # _ihm_model_group_link
            new_lines.append("loop_\n")
            new_lines.append("_ihm_model_group_link.group_id\n"
                             "_ihm_model_group_link.model_id\n")
            for entry in model_numbers:
                new_lines.append("1 %s\n" % (entry))
            new_lines.append("#\n#\n")
            # add the atom_site entries
            new_lines.append("loop_\n")
            # add the header
            for entry in atom_site_header_entries:
                new_lines.append("%s\n" % (entry))
            # and the entries
            for entry in atom_site_entries:
                new_lines.append("%s" % (entry))
            new_lines.append("#\n")

            # write the new file
            new_atom_site_filename = \
                "temp_%s_mod.cif" % (atom_site_filename.split(".cif")[0])
            outfile = open(new_atom_site_filename, 'w')
            for entry in new_lines:
                outfile.write(entry)
            outfile.close()
            print("\t>>> Generated temporary file \'%s\'" % (
                new_atom_site_filename))
            return new_atom_site_filename

def main(excel_filename, cifout_filename, atom_site_filename):
    # For debugging: use default values
    if DEBUG:
        excel_filename = 'excel_example.xlsx'
        cifout_filename = 'sample_output.cif'
        atom_site_filename = 'atom_site_input_example.cif'
    this_excel2mmcif_converter = Excel2mmcifConverter(
        excel_filename=excel_filename,
        cifout_filename=cifout_filename,
        atom_site_filename=atom_site_filename)
    this_excel2mmcif_converter.run()

if __name__ == '__main__':
    # For debugging, make the parameters not required
    parameter_requirement = True
    excel_filename = None
    cifout_filename = None
    atom_site_filename = None
    if DEBUG:
        parameter_requirement = False
    ## Command line parameters
    parser = argparse.ArgumentParser(
        description='Read an excel file and convert it the information '
                    'to mmcif.')
    parser.add_argument(
        '--excel', '-e', nargs=1, help='Name of the excel file to be read.',
        required=parameter_requirement)
    parser.add_argument(
        '--cif', '-c', nargs=1, help='Name of the mmcif file to be generated.',
        required=parameter_requirement)
    parser.add_argument(
        '--atoms', '-a', nargs=1, required=parameter_requirement,
        help='Name of the file containing the _atom_site entries that will '
             'be read. Should be in mmcif format.')
    args = parser.parse_args()

    if not DEBUG:
        excel_filename = args.excel[0]
        cifout_filename = args.cif[0]
        atom_site_filename = args.atoms[0]

    main(excel_filename, cifout_filename, atom_site_filename)
