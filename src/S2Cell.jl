
module S2Cell
#
# s2cell exceptions
#
using Markdown

md"""Exception type for invalid cell IDs."""
struct InvalidCellID <: Exception end


md"""Exception type for invalid tokens."""
struct InvalidToken <: Exception end

md"""An S2 cell ID is an integer we package for type-checing"""
struct Cell
    id::UInt32
end

#
# S2 base constants needed for cell mapping
#

# The maximum level supported within an S2 cell ID. Each level is represented by two bits in the
# final cell ID
_S2_MAX_LEVEL = 30

# The maximum value within the I and J bits of an S2 cell ID
_S2_MAX_SIZE = 1 << _S2_MAX_LEVEL

# The number of bits in a S2 cell ID used for specifying the base face
_S2_FACE_BITS = 3

# The number of bits in a S2 cell ID used for specifying the position along the Hilbert curve
_S2_POS_BITS = 2 * _S2_MAX_LEVEL + 1

# The maximum value of the Si/Ti integers used when mapping from IJ to ST. This is twice the max
# value of I and J, since Si/Ti allow referencing both the center and edge of a leaf cell
_S2_MAX_SI_TI = 1 << (_S2_MAX_LEVEL + 1)

# Mask that specifies the swap orientation bit for the Hilbert curve
_S2_SWAP_MASK = 1

# Mask that specifies the invert orientation bit for the Hilbert curve
_S2_INVERT_MASK = 2

# The number of bits per I and J in the lookup tables
_S2_LOOKUP_BITS = 4

# Lookup table for mapping 10 bits of IJ + orientation to 10 bits of Hilbert curve position +
# orientation. 
_S2_LOOKUP_POS = nothing

# Lookup table for mapping 10 bits of Hilbert curve position + orientation to 10 bits of IJ +
# orientation. 
_S2_LOOKUP_IJ = nothing

# Lookup table of two bits of IJ from two bits of curve position, based also on the current curve
# orientation from the swap and invert bits
_S2_POS_TO_IJ = [
    [0, 1, 3, 2],  # 0: Normal order, no swap or invert
    [0, 2, 3, 1],  # 1: Swap bit set, swap I and J bits
    [3, 2, 0, 1],  # 2: Invert bit set, invert bits
    [3, 1, 0, 2],  # 3: Swap and invert bits set
]

# Lookup for the orientation update mask of one of the four sub-cells within a higher level cell.
# This mask is XOR'ed with the current orientation to get the sub-cell orientation.
_S2_POS_TO_ORIENTATION_MASK = [_S2_SWAP_MASK, 0, 0, _S2_SWAP_MASK | _S2_INVERT_MASK]


#
# S2 helper functions
#

md"""
    Convert S2 UV to ST.

    This is done using the quadratic projection that is used by default for S2. The C++ and Java S2
    libraries use a different definition of the ST cell-space, but the end result in IJ is the same.
    The below uses the C++ ST definition.

    See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L317-L320

"""
function _s2_uv_to_st(component:: Real) :: Real
    if component >= 0.0
        return 0.5 * sqrt(1.0 + 3.0 * component)
    else
        return 1.0 - 0.5 * math.sqrt(1.0 - 3.0 * component)
    end
end


md"""
Convert S2 ST to UV.

This is done using the quadratic projection that is used by default for S2. The C++ and Java S2
libraries use a different definition of the ST cell-space, but the end result in IJ is the same.
The below uses the C++ ST definition.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L312-L315

"""

function _s2_st_to_uv(component:: Real) :: Real
    if component >= 0.5
        return (1.0 / 3.0) * (4.0 * component^2 - 1.0)
    else
        return (1.0 / 3.0) * (1.0 - 4.0 * (1.0 - component)^2)
    end
end

md"""
Convert S2 ST to IJ.

The mapping here differs between C++ and Java versions, but the combination of
_st_to_ij(_uv_to_st(val)) is the same for both. The below uses the C++ ST definition.

See s2geometry/blob/2c02e21040e0b82aa5719e96033d02b8ce7c0eff/src/s2/s2coords.h#L333-L336

"""
function _s2_st_to_ij(component:: Real) :: Int
    # The reference implementation does round(_S2_MAX_SIZE * component - 0.5), which is equivalent
    # to math.floor(_S2_MAX_SIZE * component)
    return int(max(0, min(_S2_MAX_SIZE - 1, math.floor(_S2_MAX_SIZE * component))))
end

md"""
Convert S2 Si/Ti to ST.

This converts an integer in range 0 to _S2_MAX_SI_TI into a float in range 0.0 to 1.0.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L338-L341

"""
function _s2_si_ti_to_st(component:: Int) :: Real
    return (1.0 / _S2_MAX_SI_TI) * component
end

md"""
Convert face + UV to S2Point XYZ.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L348-L357

Args:
    face: The S2 face for the input point.
    uv: The S2 face UV coordinates.

Returns:
    The unnormalised S2Point XYZ.

Raises:
    ValueError: If the face is not valid in range 0-5.

"""
function _s2_face_uv_to_xyz(face:: Int, uv:: Tuple{<:Real,<:Real}):: Tuple{<:Real,<:Real}
    # Face -> XYZ components -> indices with negation:
    # 0    -> ( 1,  u,  v)   -> ( /,  0,  1)
    # 1    -> (-u,  1,  v)   -> (-0,  /,  1)
    # 2    -> (-u, -v,  1)   -> (-0, -1,  /)
    # 3    -> (-1, -v, -u)   -> (-/, -1, -0) <- -1 here means -1 times the value in index 1,
    # 4    -> ( v, -1, -u)   -> ( 1, -/, -0)    not index -1
    # 5    -> ( v,  u, -1)   -> ( 1,  0, -/)
    if face == 0
        s2_point = (1, uv[0], uv[1])
    elseif face == 1
        s2_point = (-uv[0], 1, uv[1])
    elseif face == 2
        s2_point = (-uv[0], -uv[1], 1)
    elseif face == 3
        s2_point = (-1, -uv[1], -uv[0])
    elseif face == 4
        s2_point = (uv[1], -1, -uv[0])
    elseif face == 5
        s2_point = (uv[1], uv[0], -1)
    else
        throw(ValueError("Cannot convert UV to XYZ with invalid face: {}".format(face)))
    end

    return s2_point
end


#
# Cell ID <-> Token translation functions
#

md"""
Convert S2 cell ID to a S2 token.

Converts the S2 cell ID to hex and strips any trailing zeros. The 0 cell ID token is represented
as 'X' to prevent it being an empty string.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L204-L220

Args:
    cell_id: The S2 cell ID integer.

Returns:
    The S2 token string for the S2 cell ID.

Raises:
    TypeError: If the cell_id is not int.

"""
function cell_id_to_token(cell_id:: Cell)::String
    # The zero token is encoded as 'X' rather than as a zero-length string
    if cell_id == 0
        return 'X'
    end

    # Convert cell ID to 16 character hex string and strip any implicit trailing zeros
    return "{:016x}".format(cell_id).rstrip("0")
end

md"""
Convert S2 token to S2 cell ID.

Restores the stripped 0 characters from the token and converts the hex string to integer.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L222-L239

Args:
    token: The S2 token string. Can be upper or lower case hex string.

Returns:
   The S2 cell ID for the S2 token.

Raises:
    TypeError: If the token is not str.
    InvalidToken: If the token length is over 16.

"""
function token_to_cell_id(token:: String) ::Int
    # Check input
    if not isinstance(token, String)
        throw(TypeError("Cannot convert S2 token from type: {}".format(type(token))))
    end

    if len(token) > 16
        throw(InvalidToken("Cannot convert S2 token with length > 16 characters"))
    end

    # Check for the zero cell ID represented by the character 'x' or 'X' rather than as the empty
    # string
    if lowercase(token) == "x"
        return 0
    end

    # Add stripped implicit zeros to create the full 16 character hex string
    token = token + ('0' ^ (16 - length(token)))

    # Convert to cell ID by converting hex to int
    return int(token, 16)
end

#
# Encode functions
#

md"""
Convert lat/lon to a S2 cell ID.

It is expected that the lat/lon provided are normalised, with latitude in the range -90 to 90.

Args:
    lat: The latitude to convert, in degrees.
    lon: The longitude to convert, in degrees.
    level: The level of the cell ID to generate, from 0 up to 30.

Returns:
    The S2 cell ID for the lat/lon location.

Raises:
    ValueError: When level is not an integer, is < 0 or is > 30.

"""
function lat_lon_to_cell_id(lat:: Real, lon:: Real, level:: Int = 30)::Int
    if not isinstance(level, int) || level < 0 || level > _S2_MAX_LEVEL
        throw(ValueError("S2 level must be integer >= 0 and <= 30"))
    end

    # Populate _S2_LOOKUP_POS on first run.
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L75-L109
    #
    # This table takes 10 bits of I and J and orientation and returns 10 bits of curve position and
    # new orientation

    # Reuse constant expressions
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    sin_lat_rad = math.sin(lat_rad)
    cos_lat_rad = math.cos(lat_rad)
    sin_lon_rad = math.sin(lon_rad)
    cos_lon_rad = math.cos(lon_rad)

    # Convert to S2Point
    # This is effectively the unit non-geodetic ECEF vector
    s2_point = (
        cos_lat_rad * cos_lon_rad,  # X
        cos_lat_rad * sin_lon_rad,  # Y
        sin_lat_rad                 # Z
    )

    # Get cube face
    # See s2geometry/blob/2c02e21040e0b82aa5719e96033d02b8ce7c0eff/src/s2/s2coords.h#L380-L384
    #
    # The face is determined by the largest XYZ component of the S2Point vector. When the component
    # is negative, the second set of three faces is used.
    # Largest component -> face:
    # +x -> 0
    # +y -> 1
    # +z -> 2
    # -x -> 3
    # -y -> 4
    # -z -> 5
    face = maximum(abs.(s2_point))
    if s2_point[face] < 0.0
        face += 3
    end

    # Convert face + XYZ to cube-space face + UV
    # See s2geometry/blob/2c02e21040e0b82aa5719e96033d02b8ce7c0eff/src/s2/s2coords.h#L366-L372
    #
    # The faces are oriented to ensure continuity of curve.
    # Face -> UV components -> indices with negation (without divisor, which is always the remaining
    # component (index: face % 3)):
    # 0 -> ( y,  z) -> ( 1,  2)
    # 1 -> (-x,  z) -> (-0,  2)
    # 2 -> (-x, -y) -> (-0, -1) <- -1 here means -1 times the value in index 1, not index -1
    # 3 -> ( z,  y) -> ( 2,  1)
    # 4 -> ( z, -x) -> ( 2, -0)
    # 5 -> (-y, -x) -> (-1, -0)
    #
    # For a compiled language, a switch statement on face is preferable as it will be more easily
    # optimised as a jump table etc; but in Python the indexing method is more concise.
    #
    # The index selection can be reduced to some bit magic:
    # U: 1 - ((face + 1) >> 1)
    # V: 2 - (face >> 1)
    #
    # The negation of the the two components is then selected:
    # U: (face in [1, 2, 5]) ? -1: 1
    # V: (face in [2, 4, 5])) ? -1: 1
    uv = (
        s2_point[1 - ((face + 1) >> 1)] / s2_point[face % 3],  # U
        s2_point[2 - (face >> 1)] / s2_point[face % 3]         # V
    )
    if face in (1, 2, 5)
        uv = (-uv[0], uv[1])  # Negate U
    end
    if face in (2, 4, 5) 
        uv = (uv[0], -uv[1])  # Negate V
    end

    # Project cube-space UV to cell-space ST
    # See s2geometry/blob/2c02e21040e0b82aa5719e96033d02b8ce7c0eff/src/s2/s2coords.h#L317-L320
    st = (_s2_uv_to_st(uv[0]), _s2_uv_to_st(uv[1]))  # pylint: disable=invalid-name

    # Convert ST to IJ integers
    # See s2geometry/blob/2c02e21040e0b82aa5719e96033d02b8ce7c0eff/src/s2/s2coords.h#L333-L336
    ij = (_s2_st_to_ij(st[0]), _s2_st_to_ij(st[1]))  # pylint: disable=invalid-name

    # Convert face + IJ to cell ID
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L256-L298
    #
    # This is done by looking up 8 bits of I and J (4 each) at a time in the lookup table, along
    # with two bits of orientation (swap (1) and invert (2)). This gives back 8 bits of position
    # along the curve and two new orientation bits for the curve within the sub-cells in the next
    # step.
    #
    # The swap bit swaps I and J with each other
    # The invert bit inverts the bits of I and J, which means axes are negated
    #
    # Compared to the standard versions, we check the required number of steps we need to do for the
    # requested level and don't perform steps that will be completely overwritten in the truncation
    # below, rather than always doing every step. Each step does 4 bits each of I and J, which is 4
    # levels, so the required number of steps is ceil((level + 2) / 4), when level is > 0. The
    # additional 2 levels added are required to account for the top 3 bits (4 before right shift)
    # that are occupied by the face bits.
    bits = face & _S2_SWAP_MASK  # iiiijjjjoo. Initially set by by face
    cell_id = face << (_S2_POS_BITS - 1)  # Insert face at second most signficant bits
    lookup_mask = (1 << int(_S2_LOOKUP_BITS)) - 1  # Mask of 4 one bits: 0b1111
    if level > 0
        required_steps = math.ceil((level + 2) / 4)
    else
        required_steps = 0
    end
    for k in 7:-1:(7 - required_steps)
        # Grab 4 bits of each of I and J
        offset = k * _S2_LOOKUP_BITS
        bits += ((ij[0] >> offset) & lookup_mask) << (_S2_LOOKUP_BITS + 2)
        bits += ((ij[1] >> offset) & lookup_mask) << 2

        # Map bits from iiiijjjjoo to ppppppppoo using lookup table
        bits = _S2_LOOKUP_POS[bits]

        # Insert position bits into cell ID
        cell_id |= (bits >> 2) << (k * 2 * _S2_LOOKUP_BITS)

        # Remove position bits, leaving just new swap and invert bits for the next round
        bits &= _S2_SWAP_MASK | _S2_INVERT_MASK  # Mask: 0b11
    end
        
    # Left shift and add trailing bit
    # The trailing bit addition is disabled, as we are overwriting this below in the truncation
    # anyway. This line is kept as an example of the full method for S2 cell ID creation as is done
    # in the standard library versions.
    cell_id = cell_id << 1  # + 1

    # Truncate to desired level
    # This is done by finding the mask of the trailing 1 bit for the specified level, then zeroing
    # out all bits less significant than this, then finally setting the trailing 1 bit. This is
    # still necessary to do even after a reduced number of steps `required_steps` above, since each
    # step contains multiple levels that may need partial overwrite. Additionally, we need to add
    # the trailing 1 bit, which is not yet set above.
    least_significant_bit_mask = 1 << (2 * (_S2_MAX_LEVEL - level))
    cell_id = (cell_id & -least_significant_bit_mask) | least_significant_bit_mask

    return cell_id
end

md"""
Convert lat/lon to a S2 token.

Converts the S2 cell ID to hex and strips any trailing zeros. The 0 cell ID token is represented
as 'X' to prevent it being an empty string.

It is expected that the lat/lon provided are normalised, with latitude in the range -90 to 90.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L204-L220

Args:
    lat: The latitude to convert, in degrees.
    lon: The longitude to convert, in degrees.
    level: The level of the cell ID to generate, from 0 up to 30.

Returns:
    The S2 token string for the lat/lon location.

Raises:
    ValueError: When level is not an integer, is < 0 or is > 30.

"""
function lat_lon_to_token(lat:: Real, lon:: Real, level:: Int = 30)::String
    # Generate cell ID and convert to token
    return cell_id_to_token(lat_lon_to_cell_id(lat, lon, level))
end

#
# Decode functions
#

md"""
Convert S2 cell ID to lat/lon.

Args:
    cell_id: The S2 cell ID integer.

Returns:
    The lat/lon (in degrees) tuple generated from the S2 cell ID.

Raises:
    TypeError: If the cell_id is not int.
    InvalidCellID: If the cell_id is invalid.

"""
function cell_id_to_lat_lon(cell_id:: Int):: Tuple
    # Check input
    if not cell_id_is_valid(cell_id)
        throw(InvalidCellID("Cannot decode invalid S2 cell ID: {}".format(cell_id)))
    end

    # Populate _S2_LOOKUP_IJ on first run.
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L75-L109
    # This table takes 10 bits of curve position and orientation and returns 10 bits of I and J and
    # new orientation

    # Extract face + IJ from cell ID
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L312-L367
    #
    # This is done by looking up 8 bits of curve position at a time in the lookup table, along with
    # two bits of orientation (swap (1) and invert (2)). This gives back 8 bits of I and J (4 each)
    # and two new orientation bits for the curve within the sub-cells in the next step.
    #
    # The swap bit swaps I and J with each other
    # The invert bit inverts the bits of I and J, which means axes are negated
    #
    # In the first loop (most significant bits), the 3 bits occupied by the face need to be masked
    # out, since these are not set in the IJ to cell ID during encoding.
    #
    # The I and J returned here are of one of the two leaf (level 30) cells that are located
    # diagonally closest to the cell center. This happens because repeated ..00.. will select the
    # 'lower left' (for nominally oriented Hilbert curve segments) of the sub-cells. The ..10..
    # arising from the trailing bit, prior to the repeated ..00.. bits, ensures we first pick the
    # 'upper right' of the cell, then iterate in to lower left until we hit the leaf cell. This
    # means we pick the leaf cell to the north east of the parent cell center (again for nominal
    # orientation).
    # However, in the case of the swapped and inverted curve segment (4th sub-curve segment), the
    # ..10.. will select the 'lower left' and then iterate to the 'upper right' with each ..00..
    # following. In that case, we will be offset left and down by one leaf cell in each of I and J,
    # which needs to be fixed to have a consistent mapping. This is detectable by seeing that the
    # final bit of I or J is 1 (i.e we have picked an odd row/column, which will happen concurrently
    # in both I and J, so we only need to check one), except in case of level 29 where the logic is
    # inverted and the correction needs to be applied when we pick an even row/column (i.e I/J ends
    # in 0), since there are no trailing ..00..  available after the ``..10..`` when we are at level
    # 29+.
    #
    # This behaviour can be captured in the expression:
    # apply_correction = not leaf and (i ^ (is level 29)) & 1
    # apply_correction = not leaf and (i ^ (cell_id >> 2)) & 1
    #
    # We check for level 29 by looking for the trailing 1 in the third LSB, when we already know
    # that we are not a leaf cell (which could give false positive) by the initial check in the
    # expression.
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.h#L503-L529
    #
    face = cell_id >> _S2_POS_BITS
    bits = face & _S2_SWAP_MASK  # ppppppppoo. Initially set by by face
    lookup_mask = (1 << _S2_LOOKUP_BITS) - 1  # Mask of 4 one bits: 0b1111
    i = 0
    j = 0
    for k in 7:-1:-1
        # Pull out 8 bits of cell ID, except in first loop where we pull out only 4
        if k == 7
            n_bits = (_S2_MAX_LEVEL - 7 * _S2_LOOKUP_BITS)
        else
            n_bits = _S2_LOOKUP_BITS
        end
        extract_mask = (1 << (2 * n_bits)) - 1  # 8 (or 4) one bits
        bits += ((cell_id >> (k * 2 * _S2_LOOKUP_BITS + 1)) & extract_mask) << 2

        # Map bits from ppppppppoo to iiiijjjjoo using lookup table
        bits = _S2_LOOKUP_IJ[bits]

        # Extract I and J bits
        offset = k * _S2_LOOKUP_BITS
        i += (bits >> (_S2_LOOKUP_BITS + 2)) << offset  # Don't need lookup mask here
        j += ((bits >> 2) & lookup_mask) << offset

        # Remove I and J bits, leaving just new swap and invert bits for the next round
        bits &= _S2_SWAP_MASK | _S2_INVERT_MASK  # Mask: 0b11
    end

    # Resolve the center of the cell. For leaf cells, we add half the leaf cell size. For non-leaf
    # cells, we currently have one of either two cells diagonally around the cell center and want
    # to pick the leaf-cell edges that represent the parent cell center, as described above. The
    # center_correction_delta is 2x the offset, as we left shift I and J first.
    # This gives us the values Si and Ti, which are discrete representation of S and T in range 0 to
    # _S2_MAX_SI_TI. The extra power of 2 over IJ allows for identifying both the center and edge of
    # cells, whilst IJ is just the leaf cells.
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L57-L65
    is_leaf = (cell_id & 1)  # Cell is leaf cell when trailing one bit is in LSB
    apply_correction = !is_leaf & ((i ^ (cell_id >> 2)) & 1)
    if is_leaf
        correction_delta = 1
    else
        correction_delta = apply_correction ? 2 : 0
    end
    si = (i << 1) + correction_delta  # pylint: disable=invalid-name
    ti = (j << 1) + correction_delta  # pylint: disable=invalid-name

    # Convert integer si/ti to double ST
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L338-L341
    st = (_s2_si_ti_to_st(si), _s2_si_ti_to_st(ti))  # pylint: disable=invalid-name

    # Project cell-space ST to cube-space UV
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L312-L315
    uv = (_s2_st_to_uv(st[0]), _s2_st_to_uv(st[1]))  # pylint: disable=invalid-name

    # Convert face + UV to S2Point XYZ

    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2coords.h#L348-L357
    s2_point = _s2_face_uv_to_xyz(face, uv)

    # Normalise XYZ S2Point vector
    # This section is part of the reference implementation but is not necessary when mapping
    # straight into lat/lon, since the normalised and unnormalised triangles used to calculate the
    # angles are geometrically similar. If anything, the normalisation process loses precision when
    # tested against the reference implementation, albeit not at a level that is important either
    # way. The code below is left for demonstration of the normalisation process.
    # norm = math.sqrt(s2_point[0]^2 + s2_point[1]^2 + s2_point[2]^2)
    # s2_point = (s2_point[0] / norm, s2_point[1] / norm, s2_point[2] / norm)

    # Map into lat/lon
    # See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2latlng.h#L196-L205
    lat_rad = math.atan2(s2_point[2], sqrt(s2_point[0]^2 + s2_point[1]^2))
    lon_rad = math.atan2(s2_point[1], s2_point[0])

    return (math.degrees(lat_rad), math.degrees(lon_rad))
end

md"""
Convert S2 token to lat/lon.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L222-L239

Args:
    token: The S2 token string. Can be upper or lower case hex string.

Returns:
    The lat/lon (in degrees) tuple generated from the S2 token.

Raises:
    TypeError: If the token is not str.
    InvalidToken: If the token length is over 16.
    InvalidToken: If the token is invalid.
    InvalidCellID: If the contained cell_id is invalid.

"""
function token_to_lat_lon(token:: String)::Tuple
    # Check input
    if !token_is_valid(token)
        throw(InvalidToken("Cannot decode invalid S2 token: {}".format(token)))
    end

    # Convert to cell ID and decode to lat/lon
    return cell_id_to_lat_lon(token_to_cell_id(token))
end

#
# Token canonicalisation
#

md"""
Convert S2 token to a canonicalised S2 token.

This produces a token that matches the form generated by the reference C++ implementation:

- Lower case (except 'X' below)
- No whitespace
- Trailing '0' characters stripped
- Zero cell ID represented as 'X', not 'x' or ''

Args:
    token: The S2 token string to canonicalise.

Returns:
    The canonicalised S2 token.

"""
function token_to_canonical_token(token:: String):: String
    # Convert token to lower case.
    # Note that 'X' below will be returned upper case
    token = token.lower()

    # Strip any surrounding whitespace
    token = token.strip()

    # Strip any trailing zeros
    token = token.rstrip("0")

    # If empty string or 'x', return 'X' token
    if token in ("", "x")
        token = 'X'
    end

    return token
end

#
# Validation
#

md"""
Check that a S2 cell ID is valid.

Looks for valid face bits and a trailing 1 bit in one of the correct locations.

Args:
    cell_id: The S2 cell integer to validate.

Returns:
    True if the cell ID is valid, False otherwise.

Raises:
    TypeError: If the cell_id is not int.

"""
function cell_id_is_valid(cell_id:: Int):: Bool
    # Check input
    if ! cell_id isa Int
        throw(TypeError("Cannot decode S2 cell ID from type: {}".format(type(cell_id))))
    end

    # Check for zero ID
    # This avoids overflow warnings below when 1 gets added to max uint64
    if cell_id == 0
        return false
    end

    # Check face bits
    if (cell_id >> _S2_POS_BITS) > 5
        return false
    end

    # Check trailing 1 bit is in one of the even bit positions allowed for the 30 levels, using the
    # mask: 0b0001010101010101010101010101010101010101010101010101010101010101 = 0x1555555555555555
    lowest_set_bit = cell_id & (¬cell_id + 1)  # pylint: disable=invalid-unary-operand-type
    if not lowest_set_bit & 0x1555555555555555
        return false
    end

    return true  # Checks have passed, cell ID must be valid
end


md"""
Check that a S2 token is valid.

Looks for valid characters, then checks that the contained S2 cell ID is also valid. Note that
the '', 'x' and 'X' tokens are considered invalid, since the cell IDs they represent are
invalid.

Args:
    token: The S2 token string to validate.

Returns:
    True if the token is valid, False otherwise.

Raises:
    TypeError: If the token is not str.

"""
function token_is_valid(token:: String):: Bool
    # Check input
    if is(token, String)
        throw(TypeError("Cannot check S2 token with type: {}".format(type(token))))
    end

    # First check string with regex
    if !re.match(r"^[0-9a-fA-f]{1,16}$", token)
        return false
    end

    # Check the contained cell ID is valid
    return cell_id_is_valid(token_to_cell_id(token))
end

#
# Level extraction functions
#

md"""
Get the level for a S2 cell ID.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.h#L543-L551

Args:
    cell_id: The S2 cell ID integer.

Returns:
    The level of the S2 cell ID.

Raises:
    TypeError: If the cell_id is not int.
    InvalidCellID: If the cell_id is invalid.

"""
function cell_id_to_level(cell_id:: Int):: Int
    # Check input
    if ! cell_id_is_valid(cell_id)
        throw(InvalidCellID("Cannot decode invalid S2 cell ID: {}".format(cell_id)))
    end

    # Find the position of the lowest set one bit, which will be the trailing one bit. The level is
    # given by the max level (30) minus the floored division by two of the position of the lowest
    # set bit.
    #
    # The position of the lowest set bit is found using 'count trailing zeros', which would be
    # equivalent to the C++20 function std::countr_zero() or the ctz instruction.
    lsb_pos = 0
    while cell_id != 0
        if cell_id & 1
            break
        end
        lsb_pos += 1
        cell_id >>= 1
    end

    return int(_S2_MAX_LEVEL - (lsb_pos >> 1))
end

md"""
Get the level for a S2 token.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.h#L543-L551

Args:
    token: The S2 token string. Can be upper or lower case hex string.

Returns:
    The level of the S2 token.

Raises:
    TypeError: If the token is not str.
    InvalidToken: If the token length is over 16.
    InvalidToken: If the token is invalid.
    InvalidCellID: If the contained cell_id is invalid.

"""
function token_to_level(token:: String):: Int
    # Check input
    if !token_is_valid(token)
        throw(InvalidToken("Cannot decode invalid S2 token: {}".format(token)))
    end

    # Convert to cell ID and get the level for that
    return cell_id_to_level(token_to_cell_id(token))
end

#
# Parent functions
#

md"""
Get the parent cell ID of a S2 cell ID.

Args:
    cell_id: The S2 cell ID integer.
    level: The parent level to get the cell ID for. Must be less than or equal to the current
        level of the provided cell ID. If unspecified, or None, the direct parent cell ID will
        be returned.

Returns:
    The parent cell ID at the specified level.

Raises:
    TypeError: If the cell_id is not int.
    InvalidCellID: If the cell_id is invalid.
    ValueError: If cell ID is already level 0 and level is None.
    ValueError: When level is not an integer, is < 0 or is > 30.
    ValueError: If level is greater than the provided cell ID level.

"""
function cell_id_to_parent_cell_id(cell_id:: Int, level::Int = -1)::Int
    # Check input
    if ! cell_id_is_valid(cell_id)
        throw(InvalidCellID("Cannot decode invalid S2 cell ID: {}".format(cell_id)))
    end

    # Get current level of the cell ID and check it is suitable with the requested level
    current_level = cell_id_to_level(cell_id)
    if level == -1 && current_level == 0
        throw(ValueError("Cannot get parent cell ID of a level 0 cell ID"))
    end
    if level == -1
        level = current_level - 1
    end

    if !isa(level, Int) || level < 0 || level > _S2_MAX_LEVEL
        throw(ValueError("S2 level must be integer >= 0 and <= 30"))
    end

    if level > current_level
        throw(ValueError("Cannot get level {} parent cell ID of cell ID with level {}".format(
            level, current_level
        )))
    end
    if level == current_level
        # Requested parent level is current level, return cell ID itself
        return cell_id
    end

    # Truncate to desired level
    # This is done by finding the mask of the trailing 1 bit for the specified level, then zeroing
    # out all bits less significant than this, then finally setting the trailing 1 bit. This is
    # still necessary to do even after a reduced number of steps `required_steps` above, since each
    # step contains multiple levels that may need partial overwrite. Additionally, we need to add
    # the trailing 1 bit, which is not yet set above.
    least_significant_bit_mask = 1 << (2 * (_S2_MAX_LEVEL - level))
    cell_id = (cell_id & -least_significant_bit_mask) | least_significant_bit_mask

    return cell_id
end

md"""
Get the parent token of a S2 token.

Args:
    token: The S2 token string. Can be upper or lower case hex string.
    level: The parent level to get the token for. Must be less than or equal to the current
        level of the provided toke. If unspecified, or None, the direct parent token will be
        returned.

Returns:
    The parent token at the specified level.

Raises:
    TypeError: If the token is not str.
    InvalidToken: If the token length is over 16.
    InvalidToken: If the token is invalid.
    InvalidCellID: If the contained cell_id is invalid.
    ValueError: If token is already level 0 and level is None.
    ValueError: When level is not an integer, is < 0 or is > 30.
    ValueError: If level is greater than the provided token level.

"""
function token_to_parent_token(token:: String, level:: Int = nothing)::String
    # Check input
    if ! token_is_valid(token)
        throw(InvalidToken("Cannot decode invalid S2 token: {}".format(token)))
    end
    
    # Convert to cell ID and get parent and convert back to token
    return cell_id_to_token(cell_id_to_parent_cell_id(token_to_cell_id(token), level))
end

md"""
Initialise the S2 lookups in global vars _S2_LOOKUP_POS and _S2_LOOKUP_IJ.

This generates 4 variations of a 4 level deep Hilbert curve, one for each swap/invert bit
combination. This allows mapping between 8 bits (+2 orientation) of Hilbert curve position and 8
bits (+2 orientation) of I and J, and vice versa. The new orientation bits read from the mapping
tell us the base orientation of the curve segments within the next deeper level of sub-cells.

This implementation differs in structure from the reference implementation, since it is
iterative rather than recursive. The end result is the same lookup table.

See s2geometry/blob/c59d0ca01ae3976db7f8abdc83fcc871a3a95186/src/s2/s2cell_id.cc#L75-L109

"""
function __init__()
    #global _S2_LOOKUP_POS, _S2_LOOKUP_IJ  # pylint: disable=global-statement
    # Initialise empty lookup tables
    lookup_length = 1 << (2 * _S2_LOOKUP_BITS + 2)  # = 1024
    _S2_LOOKUP_POS = [0 for i in 1:lookup_length]
    _S2_LOOKUP_IJ = [0 for i in 1:lookup_length]

    # Generate lookups for each of the base orientations given by the swap and invert bits
    # 0-3 effectively
    for base_orientation in [0, _S2_SWAP_MASK, _S2_INVERT_MASK, _S2_SWAP_MASK | _S2_INVERT_MASK]
        # Walk the 256 possible positions within a level 4 curve. This implementation is not
        # the fastest since it does not reuse the common ancestor of neighbouring positions, but
        # is simpler to read
        for pos in 0:(4^4 - 1)  # 4 levels of sub-divisions
            ij = 0  # Has pattern iiiijjjj, not ijijijij 
            orientation = base_orientation

            # Walk the pairs of bits of pos, from most significant to least, getting IJ and
            # orientation as we go
            for bit_pair_offset in 0:3
                # Bit pair is effectively the sub-cell index
                bit_pair = (pos >> ((3 - bit_pair_offset) * 2)) & 0b11
                
                # Get the I and J for the sub-cell index. These need to be spread into iiiijjjj
                # by inserting as bit positions 4 and 0
                ij_bits = _S2_POS_TO_IJ[orientation + 1][bit_pair + 1]
                ij = ( 
                    (ij << 1)  # Free up position 4 and 0 from old IJ
                    | ((ij_bits & 2) << 3)  # I bit in position 4
                    | (ij_bits & 1)  # J bit in position 0
                )

                # Update the orientation with the new sub-cell orientation
                orientation = orientation ^ _S2_POS_TO_ORIENTATION_MASK[bit_pair + 1]
            end
            println(orientation)
                    
            # Shift IJ and position to allow orientation bits in LSBs of lookup
            ij <<= 2
            pos <<= 2

            # Write lookups
            _S2_LOOKUP_POS[(ij | base_orientation) + 1] = pos | orientation
            _S2_LOOKUP_IJ[(pos | base_orientation) + 1] = ij | orientation
        end
    end
end

end
