import numpy as np
import struct

def lerp(t, a, b):
    return (1 - t) * a + t * b

class VoxelReader:
    def __init__(self, fname):
        with open(fname, 'rb') as f:
            self.subdivs = np.array(struct.unpack('=QQQQ', f.read(32)), dtype=np.int32)
            self.boundMin = np.array(struct.unpack('=fff', f.read(12)))
            self.boundMax = np.array(struct.unpack('=fff', f.read(12)))
            n = self.subdivs[0] * self.subdivs[1] * self.subdivs[2] * self.subdivs[3]
            print(f"Subdivs = {self.subdivs}")
            print(f"Bound min = {self.boundMin}")
            print(f"Bound max = {self.boundMax}")
            self.data = np.frombuffer(f.read(4 * n), dtype=np.float32).reshape(self.subdivs)

    def indexToCoord(self, indexX, indexY, indexZ):
        facX = (indexX + 0.5) / self.subdivs[0]
        facY = (indexY + 0.5) / self.subdivs[1]
        facZ = (indexZ + 0.5) / self.subdivs[2]
        coordX = lerp(facX, self.boundMin[0], self.boundMax[0])
        coordY = lerp(facY, self.boundMin[1], self.boundMax[1])
        coordZ = lerp(facZ, self.boundMin[2], self.boundMax[2])
        return coordX, coordY, coordZ

    def coordToIndex(self, coordX, coordY, coordZ):
        facX = (coordX - self.boundMin[0]) / (self.boundMax[0] - self.boundMin[0])
        facY = (coordY - self.boundMin[1]) / (self.boundMax[1] - self.boundMin[1])
        facZ = (coordZ - self.boundMin[2]) / (self.boundMax[2] - self.boundMin[2])
        indexX = max(0, min(int(np.round(facX * self.subdivs[0] - 0.5)), self.subdivs[0] - 1))
        indexY = max(0, min(int(np.round(facY * self.subdivs[1] - 0.5)), self.subdivs[1] - 1))
        indexZ = max(0, min(int(np.round(facZ * self.subdivs[2] - 0.5)), self.subdivs[2] - 1))
        return indexX, indexY, indexZ

    @property
    def voxelSize(self):
        return (self.boundMax - self.boundMin) / self.subdivs[0:3]
