import pandas as pd
import numpy as np
import gzip as gz

aa_triplet_to_singlet = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    'UNK': "X",
}

extra_singlet = { 'G':'1', 'A':'2', 'C':'3', 'U':'4', 'DG':'5', 'DA':'6', 'DT':'7', 'DC':'8', 'DU':'9'}
all_singlet = {}
all_singlet.update(aa_triplet_to_singlet)
all_singlet.update(extra_singlet)

class PDB:
    def __init__(self, read_HETATM=False, READ_MODELS=False):
        self.DF = None  # Placeholder for database connection or data structure
        self.read_HETATM = read_HETATM
        self.READ_MODELS = READ_MODELS
        self.df_dict = {
            'record_name':[],
            'atom_number':[],
            'atom_name':[],
            'alt_loc':[],
            'residue_name':[],
            'chain_id':[],
            'residue_number':[],
            'insertion_code':[],
            'x_coord':[],
            'y_coord':[],
            'z_coord':[],
            'occupancy':[],
            'temp_factor':[],
            'element':[],
            'charge':[]}
        
    def read(self,file_path):
        """Reads a PDB file and initializes the PDB object.
        This is a placeholder method and should be implemented based on the PDB structure.
        """
        if file_path.endswith('.cif') or file_path.endswith('.cif.gz'):
            return self.from_CIF(file_path)
        if file_path.endswith('.pdb.gz'):
            with gz.open(file_path, 'rt') as f:
                lines = f.readlines()
        elif file_path.endswith('.pdb'):
            with open(file_path, 'r') as f:
                lines = f.readlines()
        else:
            raise ValueError("Unsupported file format. Please provide a .pdb or .cif file. (zipped or unzipped)")
        # Process lines to extract relevant information
        l = '&'.join(lines)
        if not self.READ_MODELS:
            models = l.split('MODEL ')
            lines = models[0].split('&')
        
        
        atom_lines = [line.strip() for line in lines if (line.startswith("ATOM  ") or (self.read_HETATM and line.startswith("HETATM")))]
        
        column_names = ['record_name', 'atom_number', 'atom_name', 'alt_loc', 'residue_name',
                        'chain_id', 'residue_number', 'insertion_code', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'temp_factor', 'element', 'charge']
        x = {}
        for col in column_names:
            x[col] = []
        for line in atom_lines:
            x['record_name'].append(line[0:6].strip())
            x['atom_number'].append(int(line[6:11].strip()))
            x['atom_name'].append(line[12:16].strip())
            x['alt_loc'].append(line[16:17].strip())
            x['residue_name'].append(line[17:20].strip())
            x['chain_id'].append(line[21:22].strip())
            x['residue_number'].append(int(line[22:26].strip()))
            x['insertion_code'].append(line[26:27].strip())
            x['x_coord'].append(float(line[30:38].strip()))
            x['y_coord'].append(float(line[38:46].strip()))
            x['z_coord'].append(float(line[46:54].strip()))
            x['occupancy'].append(float(line[54:60].strip()))
            x['temp_factor'].append(float(line[60:66].strip()))
            x['element'].append(line[76:78].strip())
            x['charge'].append(line[78:80].strip())
        self.DF = pd.DataFrame(x)
        
        self.DF['alt_loc'] = self.DF['alt_loc'].astype(str)
        self.DF = self.DF[np.isin(self.DF['alt_loc'].values , [' ','A'])]
        return self.DF
    def from_CIF(self,file_path):
        """
        Converts a CIF object to a PDB object.
        This is a placeholder method and should be implemented based on the CIF structure.
        """
        c = CIF()
        db = c.read(file_path)
        # Change CIF column names to PDB format
        db.rename(columns={
            'group_PDB': 'record_name',
            'id': 'atom_number',
            'label_atom_id': 'atom_name',
            'label_alt_id': 'alt_loc',
            'label_comp_id': 'residue_name',
            'label_asym_id': 'chain_id',
            'label_seq_id': 'residue_number',
            'pdbx_PDB_ins_code': 'insertion_code',
            'Cartn_x': 'x_coord',
            'Cartn_y': 'y_coord',
            'Cartn_z': 'z_coord',
            'occupancy': 'occupancy',
            'B_iso_or_equiv': 'temp_factor',
            'type_symbol': 'element',
            'pdbx_formal_charge': 'charge'
        }, inplace=True)
        db.replace('.', '', inplace=True)
        db.replace('?', '', inplace=True)
        db.drop(columns=['label_entity_id', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num'], inplace=True, errors='ignore')
        self.DF = db
        return self.DF
    def write(self, output_path):
        # Write the DataFrame to a PDB file
        # following this
        '''
        field id	definition	length	format	range	string slicing (Python)
        1	"ATOM " or "HETATM"	6	{:6s}	01-06	[0:6]
        2	atom serial number	5	{:5d}	07-11	[6:11]
        3	atom name	4	{:^4s}	13-16	[12:16]
        4	alternate location indicator	1	{:1s}	17	[16:17]
        5	residue name	3	{:3s}	18-20	[17:20]
        6	chain identifier	1	{:1s}	22	[21:22]
        7	residue sequence number	4	{:4d}	23-26	[22:26]
        8	code for insertion of residues	1	{:1s}	27	[26:27]
        9	orthogonal coordinates for X (in Angstroms)	8	{:8.3f}	31-38	[30:38]
        10	orthogonal coordinates for Y (in Angstroms)	8	{:8.3f}	39-46	[38:46]
        11	orthogonal coordinates for Z (in Angstroms)	8	{:8.3f}	47-54	[46:54]
        12	occupancy	6	{:6.2f}	55-60	[54:60]
        13	temperature factor	6	{:6.2f}	61-66	[60:66]
        14	element symbol	2	{:>2s}	77-78	[76:78]
        15	charge on the atom	2	{:2s}	79-80	[78:80]
        '''
        with open(output_path, 'w') as f:
            for i in range(len(self.DF)):
                line = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(
                    self.DF.iloc[i]['record_name'],
                    self.DF.iloc[i]['atom_number'],
                    self.DF.iloc[i]['atom_name'],
                    self.DF.iloc[i]['alt_loc'],
                    self.DF.iloc[i]['residue_name'],
                    self.DF.iloc[i]['chain_id'],
                    self.DF.iloc[i]['residue_number'],
                    self.DF.iloc[i]['insertion_code'],
                    self.DF.iloc[i]['x_coord'],
                    self.DF.iloc[i]['y_coord'],
                    self.DF.iloc[i]['z_coord'],
                    self.DF.iloc[i]['occupancy'],
                    self.DF.iloc[i]['temp_factor'],
                    self.DF.iloc[i]['element'],
                    self.DF.iloc[i]['charge']
                )
                f.write(line)
        # print(f"PDB file written to {output_path}")
        return output_path
    def get_seq(self):
        """
        Returns the sequence of the PDB object.
        This is a placeholder method and should be implemented based on the PDB structure.
        """
        seq = {}
        chain_seq = []
        
        self.DF['chain_id_res'] = self.DF['chain_id'].astype(str) + '_' + self.DF['residue_number'].apply(lambda x: str(x).zfill(10)) + '_' + self.DF['residue_name'].astype(str)
        grouped = self.DF.groupby('chain_id_res')
        for name, group in grouped:
            chain_id, res_id, comp_id = name.split('_')
            if chain_id not in seq:
                seq[chain_id] = ''
            seq[chain_id] += aa_triplet_to_singlet.get(comp_id, 'X')
        return seq

class CIF:
    def __init__(self, read_HETATM=False, READ_MODELS=False):
        self.DF = None  # Placeholder for database connection or data structure
        self.read_HETATM = read_HETATM
        self.READ_MODELS = READ_MODELS
        self.df_dict = {
            'group_PDB':[], 
            'id':[], 
            'type_symbol':[], 
            'label_atom_id':[], 
            'label_alt_id':[],
            'label_comp_id':[], 
            'label_asym_id':[], 
            'label_entity_id':[], 
            'label_seq_id':[],
            'pdbx_PDB_ins_code':[], 
            'Cartn_x':[], 
            'Cartn_y':[], 
            'Cartn_z':[],
            'occupancy':[], 
            'B_iso_or_equiv':[], 
            'pdbx_formal_charge':[],
            'auth_seq_id':[], 
            'auth_comp_id':[], 
            'auth_asym_id':[], 
            'auth_atom_id':[],
            'pdbx_PDB_model_num':[]
        }
    def read(self,file_path):
        """
        Reads a CIF file and initializes the CIF object.
        This is a placeholder method and should be implemented based on the CIF structure.
        """
        if file_path.endswith('.pdb') or file_path.endswith('.pdb.gz'):
            return self.from_PDB(file_path)
        if file_path.endswith('.cif.gz'):
            with gz.open(file_path, 'rt') as f:
                lines = f.readlines()
        elif file_path.endswith('.cif'):
            with open(file_path, 'r') as f:
                lines = f.readlines()
        else:
            raise ValueError("Unsupported file format. Please provide a .cif or .cif.gz file.")
        # print(lines)
        # Process lines to extract relevant information
        column_names = []
        atom_lines = []
        atom_lines_started = False
        for line in lines:
            if line.startswith('_atom_site.'):
                if not atom_lines_started:
                    atom_lines_started = True
                column_names.append(line.strip().replace('_atom_site.',''))
                # print(column_names)
            elif atom_lines_started and (line.startswith('ATOM') or (self.read_HETATM and line.startswith('HETATM'))):
                atom_lines.append(line.strip())

        # print(atom_lines)
        # '''
        # _atom_site.group_PDB 
        # _atom_site.id 
        # _atom_site.type_symbol 
        # _atom_site.label_atom_id 
        # _atom_site.label_alt_id 
        # _atom_site.label_comp_id 
        # _atom_site.label_asym_id 
        # _atom_site.label_entity_id 
        # _atom_site.label_seq_id 
        # _atom_site.pdbx_PDB_ins_code 
        # _atom_site.Cartn_x 
        # _atom_site.Cartn_y 
        # _atom_site.Cartn_z 
        # _atom_site.occupancy 
        # _atom_site.B_iso_or_equiv 
        # _atom_site.pdbx_formal_charge 
        # _atom_site.auth_seq_id 
        # _atom_site.auth_comp_id 
        # _atom_site.auth_asym_id 
        # _atom_site.auth_atom_id 
        # _atom_site.pdbx_PDB_model_num 
        # '''
        # column_names = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
        #                 'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id',
        #                 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
        #                 'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge',
        #                 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id',
        #                 'pdbx_PDB_model_num']
        x = {}
        for col in column_names:
            x[col] = []
        for line in atom_lines:
            parts = [l for l in line.split(' ') if l!= '']
            # print(parts)
    
            for i, col in enumerate(column_names):
                x[col].append(parts[i])

        # print(len(x),len(x[col]))
        self.DF = pd.DataFrame(x)
        self.DF['Cartn_x'] = self.DF['Cartn_x'].astype(float)
        self.DF['Cartn_y'] = self.DF['Cartn_y'].astype(float)
        self.DF['Cartn_z'] = self.DF['Cartn_z'].astype(float)
        self.DF['occupancy'] = self.DF['occupancy'].astype(float)
        self.DF['B_iso_or_equiv'] = self.DF['B_iso_or_equiv'].astype(float)
        self.DF['id'] = self.DF['id'].astype(int)
        self.DF['label_seq_id'] = self.DF['label_seq_id'].astype(int)
        self.DF['auth_seq_id'] = self.DF['auth_seq_id'].astype(int)
        
        self.DF = self.DF[np.isin(self.DF['label_alt_id'].values , ['.','A'])]
        if not self.READ_MODELS:
            self.DF = self.DF[self.DF['pdbx_PDB_model_num']=='1']
        
        return self.DF
    
    def from_PDB(self,file_path):
        """
        Converts a PDB object to a CIF object.
        This is a placeholder method and should be implemented based on the PDB structure.
        """
        p = PDB()
        d = p.read(file_path)
        # Change PDB column names to CIF format
        # column_names = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
        #                 'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id',
        #                 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
        #                 'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge',
        #                 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id',
        #                 'pdbx_PDB_model_num']
        d.rename(columns={
            'record_name': 'group_PDB',
            'atom_number': 'id',
            'atom_name': 'label_atom_id',
            'alt_loc': 'label_alt_id',
            'residue_name': 'label_comp_id',
            'chain_id': 'label_asym_id',
            'residue_number': 'label_seq_id',
            'insertion_code': 'pdbx_PDB_ins_code',
            'x_coord': 'Cartn_x',
            'y_coord': 'Cartn_y',
            'z_coord': 'Cartn_z',
            'occupancy': 'occupancy',
            'temp_factor': 'B_iso_or_equiv',
            'element': 'type_symbol',
            'charge': 'pdbx_formal_charge'
        }, inplace=True)
        self.DF = d
        self.DF['label_entity_id'] = '1'  # Placeholder for label_entity_id
        self.DF['auth_seq_id'] = self.DF['label_seq_id']
        self.DF['auth_comp_id'] = self.DF['label_comp_id']
        self.DF['auth_asym_id'] = self.DF['label_asym_id']
        self.DF['auth_atom_id'] = self.DF['label_atom_id']
        self.DF['pdbx_PDB_model_num'] = '1'  #
        # Reorder columns to match CIF format
        # replace ' ' with '.'
        self.DF.replace(' ', '.', inplace=True)
        self.DF.replace('', '.', inplace=True)
        self.DF = self.DF[[
            'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
            'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id',
            'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
            'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge',
            'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id',
            'pdbx_PDB_model_num'
        ]]
        # Convert data types
        self.DF['Cartn_x'] = self.DF['Cartn_x'].astype(float)
        self.DF['Cartn_y'] = self.DF['Cartn_y'].astype(float)
        self.DF['Cartn_z'] = self.DF['Cartn_z'].astype(float)
        self.DF['occupancy'] = self.DF['occupancy'].astype(float)
        self.DF['B_iso_or_equiv'] = self.DF['B_iso_or_equiv'].astype(float)
        self.DF['id'] = self.DF['id'].astype(int)
        self.DF['label_seq_id'] = self.DF['label_seq_id'].astype(int)
        self.DF['auth_seq_id'] = self.DF['auth_seq_id'].astype(int)
        self.DF['label_entity_id'] = self.DF['label_entity_id'].astype(str)
        self.DF['auth_comp_id'] = self.DF['auth_comp_id'].astype(str)
        self.DF['auth_asym_id'] = self.DF['auth_asym_id'].astype(str)
        self.DF['auth_atom_id'] = self.DF['auth_atom_id'].astype(str)
        self.DF['pdbx_PDB_model_num'] = self.DF['pdbx_PDB_model_num'].astype(str)

        return self.DF
    def write(self, output_path):
        """
        Converts the CIF object to a CIF file.
        This is a placeholder method and should be implemented based on the CIF structure.
        """
        with open(output_path, 'w') as f:
            f.write("loop_\n")
            for col in self.DF.columns:
                f.write(f"_atom_site.{col}\n")
            for i in range(len(self.DF)):
                line = ' '.join(str(self.DF[col][i]) for col in self.DF.columns)
                f.write(f"{line}\n")
    def get_seq(self):
        """
        Returns the sequence of the CIF object.
        This is a placeholder method and should be implemented based on the CIF structure.
        """
        seq = {}
        chain_seq = []
        # first = True
        self.DF['chain_id_res'] = self.DF['label_asym_id'].astype(str) + '_' + self.DF['label_seq_id'].apply(lambda x: str(x).zfill(10)) + '_' + self.DF['label_comp_id'].astype(str)
        grouped = self.DF.groupby('chain_id_res')
        for name, group in grouped:
            chain_id, res_id, comp_id = name.split('_')
            if chain_id not in seq:
                seq[chain_id] = ''
            seq[chain_id] += aa_triplet_to_singlet.get(comp_id, 'X')
        
        # for i in range(len(self.DF)):
        #     s = self.DF['label_comp_id'][i]
        #     c = self.DF['label_asym_id'][i]
        #     ind = self.DF['label_seq_id'][i]
            
        #     if s == ' ':
        #         continue
        #     if first:
        #         first = False
        #         chain_seq.append(aa_triplet_to_singlet[s])
        #         prev_chain = c
        #         prev_ind = ind
        #     if ind != prev_ind or c != prev_chain:
        #         chain_seq.append(aa_triplet_to_singlet[s])
        #         if c!=prev_chain:
        #             seq[prev_chain] = ''.join(chain_seq)
        #             chain_seq = [aa_triplet_to_singlet[s]]
        #         prev_chain = c
        #         prev_ind = ind
        # seq[prev_chain] = ''.join(chain_seq)
        return seq


# if __name__ == "__main__":
    # file_path = '/scratch/dshah/pdb_00002a2e_xyz.cif'
    # c = CIF()
    # c.read(file_path)
    # print(c.DF)

# p = PDB()
# p.read('/scratch/dshah/Density_to_seq/pdb_00001f05_xyz.pdb')
# p = PDB()
# p.read('/scratch/dshah/Density_to_seq/pdb_00001cf0_xyz.cif')
# p.to_pdb('/scratch/dshah/Density_to_seq/pdb_00001cf0_xyz.pdb')
# c = CIF()
# c.read('/scratch/dshah/Density_to_seq/pdb_00001cf0_xyz.cif')