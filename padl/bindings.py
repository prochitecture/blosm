"""The explicit Blosm FACE attribute contract, separate from PADL name lookup."""
from dataclasses import dataclass
from enum import IntEnum
import math

from .errors import PadlError


class FaceType(IntEnum):
    BUILDING = 1
    HOLE = 2
    PART = 3
    VOLUME = 4
    BUILDING_WITH_PARTS = 5
    BUILDING_ALSO_PART = 6

    @classmethod
    def is_building(cls, value):
        return value in (cls.BUILDING, cls.BUILDING_WITH_PARTS, cls.BUILDING_ALSO_PART)

    @classmethod
    def generates_volume(cls, value):
        return value in (cls.BUILDING, cls.PART, cls.BUILDING_ALSO_PART)

    @classmethod
    def scans_parts(cls, value):
        return value in (cls.BUILDING_WITH_PARTS, cls.BUILDING_ALSO_PART)

BUILDING_USE = {0: 'unknown', 1: 'apartments', 2: 'house', 3: 'office'}
ROOF_TYPE = {1: 'flat', 2: 'skillion', 3: 'gabled', 4: 'hipped'}
OPTIONAL = {'max_elevation': 'FLOAT', 'roof_height': 'FLOAT', 'levels_wall_count': 'INT',
            'base_elevation': 'FLOAT', 'wall_min_height': 'FLOAT', 'top_height': 'FLOAT',
            'level_height': 'FLOAT', 'ground_floor_height': 'FLOAT', 'ground_floor_factor': 'FLOAT'}
BRANCH_FLAGS = ('has_max_elevation', 'has_roof_height', 'has_levels_wall_count')
FACE_TYPES = {'face_type': 'INT', 'footprint_id': 'INT', 'parent_id': 'INT', 'has_holes': 'BOOLEAN',
              'padl_random_key': 'INT',
              'building_use': 'INT', 'roof_type': 'INT', 'min_elevation': 'FLOAT',
              **OPTIONAL, **{'has_' + name: 'BOOLEAN' for name in OPTIONAL}}


@dataclass(frozen=True)
class ParameterBinding:
    name: str
    data_type: str
    presence: str | None = None
    categories: dict | None = None


class FootprintBindings:
    version = 2

    def __init__(self, extra=()):
        self.parameters = {
            'building_use': ParameterBinding('building_use', 'INT', categories=BUILDING_USE),
            'roof_type': ParameterBinding('roof_type', 'INT', categories=ROOF_TYPE),
            'min_elevation': ParameterBinding('min_elevation', 'FLOAT'),
            **{k: ParameterBinding(k, v, 'has_' + k) for k, v in OPTIONAL.items()},
            **{k: ParameterBinding(k, 'BOOLEAN') for k in BRANCH_FLAGS}}
        for binding in extra:
            if binding.name in self.parameters:
                raise ValueError('Duplicate binding: ' + binding.name)
            self.parameters[binding.name] = binding

    @staticmethod
    def check(name, value, kind):
        if kind == 'BOOLEAN':
            valid = type(value) is bool
        elif kind == 'INT':
            valid = type(value) is int and -(2**31) <= value < 2**31
        else:
            try:
                valid = type(value) in (int, float) and math.isfinite(value) and abs(value) <= 3.4028234663852886e38
            except OverflowError:
                valid = False
        if not valid:
            raise PadlError(f'{name} requires a valid {kind} value; got {value!r}', code='invalid_input')

    def decode(self, row):
        """Convert one stored FACE row to a sparse logical PADL input namespace."""
        result = {}
        for name, binding in self.parameters.items():
            if binding.presence:
                flag = row.get(binding.presence, False)
                self.check(binding.presence, flag, 'BOOLEAN')
                if not flag:
                    continue
                if name not in row:
                    raise PadlError(f'{binding.presence} is true but {name} is absent', code='missing_input')
            if name == 'roof_type' and row.get(name, 0) == 0:
                continue
            if name not in row:
                if name == 'min_elevation':
                    result[name] = 0
                elif name in BRANCH_FLAGS:
                    result[name] = False
                elif name == 'building_use':
                    result[name] = 'unknown'
                continue
            value = row[name]
            self.check(name, value, binding.data_type)
            if binding.categories is not None:
                if value not in binding.categories:
                    raise PadlError(f'Unknown {name} code {value}', code='invalid_category')
                value = binding.categories[value]
            result[name] = value
        return result

    def encode(self, inputs, *, footprint_id, has_holes=False, face_type=1, parent_id=0):
        """Initialize an outer row from current supplied logical inputs only."""
        self.check('footprint_id', footprint_id, 'INT')
        if footprint_id <= 0:
            raise PadlError('Outer footprint_id must be positive', code='invalid_identity')
        self.check('has_holes', has_holes, 'BOOLEAN')
        if type(face_type) not in (int, FaceType) or not (FaceType.is_building(face_type) or face_type == FaceType.PART):
            raise PadlError('Expected a building or part face type', code='invalid_attribute')
        self.check('parent_id', parent_id, 'INT')
        if (FaceType.is_building(face_type) and parent_id != 0) or (face_type == FaceType.PART and parent_id <= 0):
            raise PadlError('Invalid footprint parent association', code='invalid_identity')
        row = {'face_type': int(face_type), 'footprint_id': footprint_id, 'parent_id': parent_id, 'has_holes': has_holes,
               'building_use': 0, 'roof_type': 0, 'min_elevation': 0,
               **{'has_' + name: False for name in OPTIONAL}}
        unknown = inputs.keys() - self.parameters.keys()
        if unknown:
            raise PadlError(f'No footprint binding for {sorted(unknown)[0]!r}', code='unknown_binding')
        for name, value in inputs.items():
            binding = self.parameters[name]
            if binding.categories is not None:
                reverse = {v: k for k, v in binding.categories.items()}
                if value not in reverse:
                    raise PadlError(f'Unsupported {name} {value!r}', code='invalid_category')
                value = reverse[value]
            self.check(name, value, binding.data_type)
            row[name] = value
            if binding.presence:
                if inputs.get(binding.presence) is False:
                    raise PadlError(f'Conflicting supplied {name} and false {binding.presence}', code='invalid_input')
                row[binding.presence] = True
        self.decode(row)
        return row

    def validate_mesh(self, mesh):
        """Validate storage before GN conversion can silently coerce its types."""
        for name, kind in FACE_TYPES.items():
            attr = mesh.attributes.get(name)
            if attr is not None and (attr.domain != 'FACE' or attr.data_type != kind):
                raise PadlError(f'{name} must be FACE {kind}', code='invalid_attribute')
        for name in ('face_type', 'footprint_id', 'parent_id', 'has_holes'):
            if mesh.attributes.get(name) is None:
                raise PadlError(f'Required FACE attribute {name} is absent', code='missing_attribute')
        rows = [{name: mesh.attributes[name].data[i].value for name in FACE_TYPES if mesh.attributes.get(name) is not None}
                for i in range(len(mesh.polygons))]
        outer = {}
        holes = {}
        for index, row in enumerate(rows):
            role = row['face_type']
            if FaceType.is_building(role) or role == FaceType.PART:
                fid = row['footprint_id']
                if fid <= 0 or fid in outer or (FaceType.is_building(role) and row['parent_id'] != 0):
                    raise PadlError('Footprint IDs must be unique and positive; buildings require parent_id=0', code='invalid_identity')
                outer[fid] = index
                self.decode(row)
            elif role == 2:
                if row['footprint_id'] != 0 or row['parent_id'] <= 0:
                    raise PadlError('Holes require footprint_id=0 and a positive parent_id', code='invalid_identity')
                holes.setdefault(row['parent_id'], []).append(index)
            elif role not in (0, 4):
                raise PadlError(f'Unknown face_type {role}', code='invalid_attribute')
        for parent in holes:
            if parent not in outer:
                raise PadlError(f'Hole has no outer footprint {parent}', code='invalid_identity')
        for fid, index in outer.items():
            row = rows[index]
            if row['face_type'] == FaceType.PART:
                parent = outer.get(row['parent_id'])
                if parent is None or not FaceType.scans_parts(rows[parent]['face_type']):
                    raise PadlError(f'Part {fid} requires a type-5 or type-6 building parent', code='invalid_identity')
            if row['face_type'] == FaceType.BUILDING_WITH_PARTS and not any(
                    r['face_type'] == FaceType.PART and r['parent_id'] == fid for r in rows):
                raise PadlError(f'Type-5 building {fid} has no associated parts', code='invalid_identity')
            if rows[index]['has_holes'] != bool(holes.get(fid)):
                raise PadlError(f'has_holes does not match the associated holes of {fid}', code='invalid_identity')
        return rows, outer, holes
