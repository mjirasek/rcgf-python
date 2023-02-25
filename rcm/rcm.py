import pandas as pd
import numpy as np
from rdkit import Chem
from typing import Iterable
from . import biot_savart
from . import area
from .messages import WarningMessage, LoadingMessage


class RCM:
    
    
    def __init__(self, 
                xyz: str, 
                conn: str) -> None:
        self.mol: Chem.rdchem.Mol = self._xyz_to_mol(xyz)
        self.xyz: pd.DateFrame = self._mol_to_pd_xyz()
        self.conn: pd.DataFrame = self._connectivity_reader(conn)
        self.M = self._combined_matrix()
        self.check_if_current_flow_conserved()

        try: 
            self.net_area: tuple[float] = self.get_net_area_weight()
        except:
            print(WarningMessage.AREA_ESTIMATION_FAIL)

    def _xyz_to_mol(self, filepath: str):
        # print(LoadingMessage.LOADING_XYZ)
        return Chem.rdmolfiles.MolFromXYZFile(filepath)

    def _combined_matrix(self) -> pd.DataFrame:
        M = pd.DataFrame(
            np.nan, 
            columns = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'J'], 
            index = np.arange(len(self.conn))
        )
        
        for i in self.conn.index:
            temp_start: int = self.conn.loc[i,'start']
            temp_end: int = self.conn.loc[i,'end']
            J: float = self.conn.loc[i,'current_weight']

            a2: float = self.xyz.loc[(temp_start-1),'x']
            b2: float = self.xyz.loc[(temp_start-1),'y']
            c2: float = self.xyz.loc[(temp_start-1),'z']

            a1: float = self.xyz.loc[(temp_end-1),'x'] - a2
            b1: float = self.xyz.loc[(temp_end-1),'y'] - b2
            c1: float = self.xyz.loc[(temp_end-1),'z'] - c2

            M.loc[i,::] = a1,a2,b1,b2,c1,c2,J
        return M

    def _connectivity_reader(self, filepath: str) -> pd.DataFrame:
        # print(LoadingMessage.LOADING_CONN)
        return (
            pd.read_csv(
            filepath, header=None, 
            names = ['current_weight', 'start','end']
            )
        )

    def _mol_to_pd_xyz(self) -> pd.DataFrame:
        def count_iterable(i: Iterable) -> int:
            return sum(1 for e in i)
        df_xyz = pd.DataFrame(
            np.nan, 
            columns=['atom', 'x', 'y', 'z'], 
            index=np.arange(count_iterable(self.mol.GetAtoms()))
        )

        for i,atom in enumerate(self.mol.GetAtoms()): 
            positions = self.mol.GetConformer().GetAtomPosition(i)
            df_xyz.at[i,'atom'] = atom.GetSymbol()
            df_xyz.at[i,'x'] = positions.x
            df_xyz.at[i,'y'] = positions.y
            df_xyz.at[i,'z'] = positions.z

        return df_xyz

    def get_B(
            self, 
            x: float, 
            y: float, 
            z: float
        ) -> tuple[np.float64]:
        
        return tuple(map(pd.DataFrame.sum,
                biot_savart.B(
                    self.M.a1, self.M.a2, 
                    self.M.b1, self.M.b2,
                    self.M.c1, self.M.c2, 
                    x, y, z, 
                    self.M.J
                )
            )
        )

##################################
#################################
    def get_np_b(self, xyz: np.array) -> np.array:
        
        return tuple(map(pd.DataFrame.sum,
                biot_savart.B(
                    self.M.a1, self.M.a2, 
                    self.M.b1, self.M.b2,
                    self.M.c1, self.M.c2, 
                    xyz[0], xyz[1], xyz[2], 
                    self.M.J
                )
            )
        )



    def _return_index_for_area_estimation(self) -> list[int]:
        no_split_path = self.conn.current_weight == 1
        if len(no_split_path) < 3:
            raise AssertionError(
                'The list of connectivities'
                ' is shorter than 3. Very'
                ' strange error.'
            )

        if sum(no_split_path) < 3:
            raise AssertionError(
                'The list of connections that'
                ' have currenth path weight 1'
                ' is shorter than 3 bonds. This'
                ' means that there are not enough'
                ' segments with current weight "1"'
                '. Thus, the area'
                ' cannot be estimated unambigiously.'
                '\n'
                'Perhaph the current paths are split'
                ' into two equal parts and both have'
                ' weight 0.5?'
            )

        def _return_row_first_bond_weight_one() -> int:
            return int(self.conn.loc[self.conn.current_weight == 1, ::].index[0])

        def _return_end_conn_number_for_row() -> int:
            return int(self.conn.loc[current_pd_index,'end'])

        def _return_current_weight() -> float:
            return float(self.conn.loc[current_pd_index,'current_weight'])

        # It starts with first segment with weight 1.
        current_pd_index = _return_row_first_bond_weight_one()
        current_conn_index = _return_end_conn_number_for_row()

        index_depo = []
        index_depo.append(current_conn_index)

        MAX_STEPS = 100
        while len(index_depo) < MAX_STEPS:
            current_pd_index = int(self.conn.loc[self.conn.start == current_conn_index, ::].index[0])
            current_conn_index = int(self.conn.loc[current_pd_index,'end'])
            if current_conn_index in index_depo:
                break
            if _return_current_weight() == 1:
                index_depo.append(current_conn_index)
        else:
            raise RuntimeWarning(
                    f'Iteration throught the connectivity '
                    f'matrix took too many steps. '
                    f'Maximum steps to walk through the '
                    f'connectivity is {MAX_STEPS}. '
                )

        return index_depo
    
    def get_net_area_weight(self) -> tuple[float]:

        # function that returns me list of indexes for which you should calculate the area, in ordered way 
        index_depo = self._return_index_for_area_estimation()
        new_index_depo = [x-1 for x in index_depo]
        self.xyz_for_area = self.xyz.loc[new_index_depo,['x','y','z']]
        Ax, Ay, Az = area.area_3d(self.xyz_for_area)
        # self.area = Ax, Ay, Az
        return Ax, Ay, Az



    @classmethod
    def grid_generator_2d(
            cls,
            x_range: float = 20, 
            y_range: float = 20,
            resolution: float = 1,
            z: float = 0,
        ) -> pd.DataFrame:
        x = np.arange(-x_range/2, x_range/2+resolution, resolution)
        y = np.arange(-y_range/2, y_range/2+resolution, resolution)
        xyz_grid = np.full([x.size*y.size, 3], np.nan)
        mesh_index = 0
        for i in x:
            for j in y:
                xyz_grid[mesh_index,::] = i,j,z
                mesh_index += 1
        xyz_grid_pd = pd.DataFrame(xyz_grid,columns=['x','y','z'])
        # xyz_grid_pd = pd.DataFrame(xyz_grid,columns=['xyz'])
        return xyz_grid_pd



    # @classmethod
    # def grid_generator_2d(
    #         cls,
    #         x_range: float = 20, 
    #         y_range: float = 20,
    #         resolution: float = 1,
    #         z: float = 0,
    #     ) -> np.array:

    #     x = np.arange(-x_range/2, x_range/2+resolution, resolution)
    #     y = np.arange(-y_range/2, y_range/2+resolution, resolution)
    #     xyz_grid = np.full([x.size*y.size, 3], np.nan)
    #     mesh_index = 0
    #     for i in x:
    #         for j in y:
    #             xyz_grid[mesh_index,::] = i,j,z
    #             mesh_index += 1
    #     return xyz_grid




    # def screen_2d(
    #         self, 
    #         x_range: float = 20, 
    #         y_range: float = 20,
    #         resolution: float = 1,
    #     ) -> None:
    #     """_summary_

    #     Args:
    #         x_range (float): 
    #                 width of the grid in Angrtom. Defaults to 20.
    #         y_range (float): 
    #                 height of the grid in Angrtom. Defaults to 20.
    #         resolution (float): 
    #                 resolution of the grid in Angstrom. Defaults to 1.
    #     """

    #     grid = self.grid_generator_2d(x_range,y_range,resolution)
    #     print(grid)
    #     # generate 2d grid
    #     # screen 2d plot
    #     pass


    def check_if_current_flow_conserved(self) -> bool:

        # first basic check if the dataframe has 3 columns.
        assert self.conn.shape[1] == 3, (
            f"The connectivity dataframe has {self.conn.shape[1]}"
            f" columns, exactly 3 are required."
        )

        # check if all columns are named as they should
        for i in ['start','end','current_weight']:
            assert i in self.conn.columns, (
            f"The connectivity dataframe does not contain" 
            f" expected column of name {i}."
            )

        # check all numbers in columns two and three
        combined_list = pd.concat([self.conn.start, self.conn.end], axis=0).unique()

        for i in combined_list:
            # check how many goes in and out
            current_in_out_balanced = (
                sum(self.conn.loc[self.conn.end == i,'current_weight']) - 
                sum(self.conn.loc[self.conn.start == i,'current_weight'])
            ) == 0
            assert current_in_out_balanced, (
                f'Current flowing in and out '
                f'of node with index {i} '
                f'is not conserved. '
            )

