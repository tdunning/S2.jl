# Copyright 2020 Adam Liddell
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using Test
using CSV
using Printf
using DataStructures

using S2Cell

@info "Low-level internals"
@test_throws ArgumentError("Cannot convert UV to XYZ with invalid face: 6") S2Cell._s2_face_uv_to_xyz(6, (0, 0))


encode_tbl = CSV.File("./s2_encode_corpus.csv.gz"; types=[Float64, Float64, Int32, UInt64, String])
decode_tbl = CSV.File("./s2_decode_corpus.csv.gz"; types=[UInt64, String, Float64, Float64, Int32])

@info "Cell to token"
@test S2Cell.cell_id_to_token(S2Cell.Cell(0)) == "X"

for row in encode_tbl
    try
        @test S2Cell.cell_id_to_token(S2Cell.Cell(UInt64(row[:cell_id]))) == row[:token]
    catch e
        println(row)
        break
    end
end

@info "Token to cell"
@test S2Cell.token_to_cell_id("x").id == 0
@test S2Cell.token_to_cell_id("X").id == 0

@test_throws S2Cell.InvalidToken("Cannot convert S2 token with length > 16 characters") S2Cell.token_to_cell_id('a' ^ 17)

for row in encode_tbl
    @test S2Cell.token_to_cell_id(row[:token]).id == row[:cell_id]
end

@info "LL to *"
@test_throws ArgumentError S2Cell.lat_lon_to_cell_id(0., 0., -1)
@test_throws ArgumentError S2Cell.lat_lon_to_cell_id(0., 0., 31)

for row in encode_tbl
    @test S2Cell.lat_lon_to_cell_id(row[:lat], row[:lon], row[:level]).id == row[:cell_id]
end

@test_throws ArgumentError S2Cell.lat_lon_to_token(0, 0, -1)
@test_throws ArgumentError S2Cell.lat_lon_to_token(0, 0, 31)

for row in encode_tbl
    @test S2Cell.lat_lon_to_token(row[:lat], row[:lon], row[:level]) == row[:token]
end

@info "* to LL"
@test_throws S2Cell.InvalidCellID S2Cell.cell_id_to_lat_lon((UInt64(0b110) << S2Cell._S2_POS_BITS))

for row in decode_tbl
    @test all(S2Cell.cell_id_to_lat_lon(row[:cell_id]) .≈ (row[:lat], row[:lon]))
end
    
@test_throws S2Cell.InvalidToken S2Cell.token_to_lat_lon('a' ^ 17)

@test_throws S2Cell.InvalidToken S2Cell.token_to_lat_lon(@sprintf "%016x" (UInt64(0b110) << S2Cell._S2_POS_BITS))

for row in decode_tbl
    @test all(S2Cell.token_to_lat_lon(row[:token]) .≈ (row[:lat], row[:lon]))
end    

@info "Token miscellany"
for (token, expected) in [
    ("3", "3"), ("2ef59bd352b93ac3", "2ef59bd352b93ac3"), ("  2ef59bd352b93ac3", "2ef59bd352b93ac3"),
    ("2ef59bd352b93ac3  ", "2ef59bd352b93ac3"), (" 2ef59bd352b93ac3 ", "2ef59bd352b93ac3"),
    ("2EF", "2ef"), ("2eF", "2ef"), ("2ef000", "2ef"), ("", "X"), ("x", "X"),
    ]
    @test S2Cell.token_to_canonical_token(token) == expected
end

evens = [ (UInt64(1) << even_number, true) for even_number in 0:2:S2Cell._S2_POS_BITS]
odds = [ (UInt64(1) << even_number, false) for even_number in 1:2:S2Cell._S2_POS_BITS]
trivial = [
    (UInt64(0b1111010101010101010101010101010101010101010101010101010101010101), false),  # Invalid face
    (UInt64(0), false),
]

for (id, is_valid) in cat(trivial, evens, odds, dims=1)
    @test S2Cell.cell_id_is_valid(id) == is_valid
end

for row in decode_tbl
    S2Cell.cell_id_is_valid(row[:cell_id])
end

for (token, is_valid) in [
    ("", false),  # Invalid cell ID
    ("x", false),  # Invalid cell ID
    ("X", false),  # Invalid cell ID
    ("2ef", true),
    ("2EF", true),
    ("2ef0000000000000", true),
    ("2ef00000000000000", false),  # Too long
    ("2efinvalid", false),  # Invalid characters
    ("2efx", false),  # Incorrect use of X
    ]
    @test S2Cell.token_is_valid(token) == is_valid
end


for row in decode_tbl
    @test S2Cell.token_is_valid(row[:token])
end

@info "Level shifting"
@test_throws S2Cell.InvalidCellID S2Cell.cell_id_to_level(0)

for row in encode_tbl
    @test S2Cell.cell_id_to_level(row[:cell_id]) == row[:level]
end

@test_throws S2Cell.InvalidToken S2Cell.token_to_level('a' ^ 17)
@test_throws S2Cell.InvalidToken S2Cell.token_to_level("")

for row in encode_tbl
    @test S2Cell.token_to_level(row[:token]) == row[:level]
end

@test_throws S2Cell.InvalidCellID S2Cell.cell_id_to_parent_cell_id(0)
@test_throws ArgumentError S2Cell.cell_id_to_parent_cell_id(3458764513820540928)
@test_throws ArgumentError S2Cell.cell_id_to_parent_cell_id(3383782026652942336, -1)
@test_throws ArgumentError S2Cell.cell_id_to_parent_cell_id(3383782026652942336, 31)
@test_throws ArgumentError S2Cell.cell_id_to_parent_cell_id(3383782026652942336, 16)

points = DefaultDict(Dict)
for row in encode_tbl
    points[(row[:lat], row[:lon])][row[:level]] = row[:cell_id]
end
for levels_dict in rand(collect(values(points)), 100)
    for level in keys(levels_dict)
        if level > 0
            @test S2Cell.cell_id_to_parent_cell_id(levels_dict[level]) == levels_dict[level-1]
        else
            @test_throws ArgumentError S2Cell.cell_id_to_parent_cell_id(levels_dict[level])
        end

        for other in keys(levels_dict)
            if other <= level
                @test S2Cell.cell_id_to_parent_cell_id(levels_dict[level], other) == levels_dict[other]
            else
                @test_throws ArgumentError S2Cell.cell_id_to_parent_cell_id(levels_dict[level], other)
            end
        end
    end
end


@test_throws S2Cell.InvalidToken S2Cell.token_to_parent_token('a' ^ 17)
@test_throws ArgumentError S2Cell.token_to_parent_token("3")
@test_throws ArgumentError S2Cell.token_to_parent_token("2ef59bd34", -1)
@test_throws ArgumentError S2Cell.token_to_parent_token("2ef59bd34", 31)
@test_throws ArgumentError S2Cell.token_to_parent_token("2ef59bd34", 16)

points = DefaultDict(Dict)
for row in encode_tbl
    points[(row[:lat], row[:lon])][row[:level]] = row[:token]
end
for levels_dict in rand(collect(values(points)), 100)
    for level in keys(levels_dict)
        if level > 0
            @test S2Cell.token_to_parent_token(levels_dict[level]) == levels_dict[level-1]
        else
            @test_throws ArgumentError S2Cell.token_to_parent_token(levels_dict[level])
        end

        for other in keys(levels_dict)
            if other <= level
                @test S2Cell.token_to_parent_token(levels_dict[level], other) == levels_dict[other]
            else
                @test_throws ArgumentError S2Cell.token_to_parent_token(levels_dict[level], other)
            end
        end
    end
end
