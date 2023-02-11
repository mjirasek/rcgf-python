import pandas as pd
import numpy as np
from rdkit import Chem
from typing import Iterable
from . import biot_savart

class RCM:
    
    
    def __init__(self, 
                xyz: str, 
                conn: str) -> None:
        self.xyz: pd.DateFrame = self._xyz_reader(xyz)
        self.conn: pd.DataFrame = self._connectivity_reader(conn)
        self.M = self._combined_matrix()
        self.check_if_current_flow_conserved()
    

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
        return (
            pd.read_csv(
            filepath, header=None, 
            names = ['current_weight', 'start','end']
            )
        )

    def _xyz_reader(self, filepath: str) -> pd.DataFrame:
        def count_iterable(i: Iterable) -> int:
            return sum(1 for e in i)
        mol = Chem.rdmolfiles.MolFromXYZFile(filepath)
        self.mol = mol
        df_xyz = pd.DataFrame(
            np.nan, 
            columns=['atom', 'x', 'y', 'z'], 
            index=np.arange(count_iterable(mol.GetAtoms()))
        )

        for i,atom in enumerate(mol.GetAtoms()): 
            positions = mol.GetConformer().GetAtomPosition(i)
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

    @classmethod
    def grid_generator_2d(
            cls,
            x_range: float = 20, 
            y_range: float = 20,
            resolution: float = 1,
            z: float = 0,
        ) -> np.array:

        x = np.arange(-x_range/2, x_range/2+resolution, resolution)
        y = np.arange(-y_range/2, y_range/2+resolution, resolution)
        xyz_grid = np.full([x.size*y.size, 3], np.nan)
        mesh_index = 0
        for i in x:
            for j in y:
                xyz_grid[mesh_index,::] = i,j,z
                mesh_index += 1
        return xyz_grid

    def screen_2d(
            self, 
            x_range: float = 20, 
            y_range: float = 20,
            resolution: float = 1,
        ) -> None:
        """_summary_

        Args:
            x_range (float): 
                    width of the grid in Angrtom. Defaults to 20.
            y_range (float): 
                    height of the grid in Angrtom. Defaults to 20.
            resolution (float): 
                    resolution of the grid in Angstrom. Defaults to 1.
        """

        grid = self.grid_generator_2d(x_range,y_range,resolution)
        print(grid)
        # generate 2d grid
        # screen 2d plot
        pass


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

