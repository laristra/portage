/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <algorithm>

#include "gtest/gtest.h"

// portage includes
#include "portage/search/search_swept_face.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

TEST(search_swept_face_2D, cell) {
  Wonton::Simple_Mesh sm{0, 0, 1, 1, 3, 3};
  Wonton::Simple_Mesh tm{0, 0, 1, 1, 3, 3};
  const Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sm);
  const Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(tm);

  Portage::SearchSweptFace<2, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_Mesh_Wrapper>
    search(source_mesh_wrapper, target_mesh_wrapper);

  // Test candidates for the interior (central) cell
  int tcellID = 4;
  std::vector<int> expected_candidates = {1, 3, 5, 7};
  std::vector<int> candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);
  
  //Test candidates for a boundary cell
  tcellID = 7;
  expected_candidates = {4, 6, 8};
  candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);

  //Test candidates for a corner cell
  tcellID = 2;
  expected_candidates = {1, 5};
  candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);
}  // TEST(search_swept_face2, cell)

TEST(search_swept_face_3D, cell) {
  Wonton::Simple_Mesh sm{0, 0, 0, 1, 1, 1, 3, 3, 3};
  Wonton::Simple_Mesh tm{0, 0, 0, 1, 1, 1, 3, 3, 3};
  const Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sm);
  const Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(tm);

  Portage::SearchSweptFace<3, Portage::Entity_kind::CELL,
                           Wonton::Simple_Mesh_Wrapper,
                           Wonton::Simple_Mesh_Wrapper>
    search(source_mesh_wrapper, target_mesh_wrapper);

  // Test candidates for the interior (central) cell
  int tcellID = 13;
  std::vector<int> expected_candidates = {4, 10, 12, 14, 16, 22};
  std::vector<int> candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);

  //Test candidates for a boundary cell
  tcellID = 16;
  expected_candidates = {7, 13, 15, 17, 25};
  candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);

  //Test candidates for another boundary cell
  tcellID = 23;
  expected_candidates = {14, 20, 22, 26};
  candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);

  //Test candidates for a corner cell
  tcellID = 0;
  expected_candidates = {1, 3, 9};
  candidates = search(tcellID);
  ASSERT_EQ(expected_candidates.size(), candidates.size());
  // candidates might not be in order, so sort them
  std::sort(candidates.begin(), candidates.end());  
  for (int icc = 0; icc < expected_candidates.size(); icc++)
    ASSERT_EQ(expected_candidates[icc], candidates[icc]);

}  // TEST(search_swept_face3, cell)