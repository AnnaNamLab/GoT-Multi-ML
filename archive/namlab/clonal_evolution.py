"""v0.1 - 2025-03-24"""


def has_vaf_conflict(path1, path2, vaf_df):
    """Check if two paths have a VAF conflict."""
    parent1 = path1[0]
    parent2 = path2[0]

    # Get the VAF values for the top parent genes
    parent1_vaf = vaf_df[vaf_df["target"] == parent1]["vaf_corrected"].values[0]
    parent2_vaf = vaf_df[vaf_df["target"] == parent2]["vaf_corrected"].values[0]

    ## Sum of top-level parent VAFs should be less than 0.5 (For case where the parents are different)
    if parent1 != parent2:
        # If sum of top-level parent VAFs is larger than 0.5, then it's a violation to VAF constraint
        if parent1_vaf + parent2_vaf > 0.5:
            return True

    ## Sum of two direct children VAFs of a common parent cannot be larger than the parent VAF
    if parent1 == parent2:
        # Loop through the children until we find different children
        for i in range(
            1, min(len(path1), len(path2))
        ):  # start from 1 because 0 is the common parent
            if path1[i] != path2[i]:
                child1 = path1[i]
                child2 = path2[i]
                child1_vaf = vaf_df[vaf_df["target"] == child1]["vaf_corrected"].values[
                    0
                ]
                child2_vaf = vaf_df[vaf_df["target"] == child2]["vaf_corrected"].values[
                    0
                ]

                # If the sum of the children VAFs is larger than 0.5, it's a violation
                if child1_vaf + child2_vaf > 0.5:
                    return True

                current_parent1 = path1[i - 1]
                current_parent2 = path2[i - 1]

                if current_parent1 == current_parent2:
                    common_parent = current_parent1
                    parent_vaf = vaf_df[vaf_df["target"] == common_parent][
                        "vaf_corrected"
                    ].values[0]

                    # If the sum of the children VAFs is larger than the parent VAF, it's a violation
                    if child1_vaf + child2_vaf > parent_vaf:
                        return True
                # if the parents (but not root) are different, the children should not have higher VAF than the parents
                elif i > 1 and current_parent1 != current_parent2:
                    current_parent1_vaf = vaf_df[vaf_df["target"] == current_parent1][
                        "vaf_corrected"
                    ].values[0]
                    current_parent2_vaf = vaf_df[vaf_df["target"] == current_parent2][
                        "vaf_corrected"
                    ].values[0]

                    if child1_vaf > current_parent1_vaf:
                        return True
                    if child2_vaf > current_parent2_vaf:
                        return True
                else:
                    raise ValueError("Invalid case")
    return False


def has_homoplasy_conflict(path1, path2):
    """Check if two paths have a homoplasy conflict (convergence evolution)."""
    # One of the paths is length 1 (One is included in the other)
    if len(path1) == 1 or len(path2) == 1:
        if path1[0] == path2[0]:  # Same top-level parent
            return False  # Allow it

        # Different top-level parent but one is included in the other, it's a homoplasy conflict
        # because the top-level parent of one path is a child of the other path
        if path1[0] != path2[0]:
            if set(path1).intersection(set(path2)):
                return True

    # Check if 2 paths have the same top-level node
    for i in range(min(len(path1), len(path2))):
        current_level_node1 = path1[i]
        current_level_node2 = path2[i]

        # If they have the same top-level node, check if they have the same next-level node
        if current_level_node1 == current_level_node2:
            next_level_node1 = path1[i + 1] if i + 1 < len(path1) else None
            next_level_node2 = path2[i + 1] if i + 1 < len(path2) else None

            # There's no next-level node for one of the paths, then it's included in the other
            if next_level_node1 is None:
                return False
            elif next_level_node2 is None:
                return False

            ## if they have the same next-level node, go to the next level
            if next_level_node1 == next_level_node2:
                continue
            ## if they have different next-level node, check if there are any common nodes from the next level (excluding the current node)
            elif next_level_node1 != next_level_node2:
                common_nodes = set(path1[i + 1 :]).intersection(set(path2[i + 1 :]))
                if common_nodes:
                    return True
                else:
                    return False

        # If they have different top-level node, check if there are any common nodes among the 2 paths
        elif current_level_node1 != current_level_node2:
            common_nodes = set(path1).intersection(set(path2))
            if common_nodes:
                return True
            else:
                return False

        else:
            raise ValueError("Invalid case")


# Function to reassign genotypes according to evolutionary trajectories
def reassign_genotypes_by_trajectory(cell_gene_matrix, path_group):
    """
    Reassign genotypes to be consistent with the evolutionary trajectories in the path group.

    Parameters:
    - cell_gene_matrix: DataFrame with cell-by-gene matrix
    - path_group: Dictionary containing path group information

    Returns:
    - reassigned_matrix: Copy of cell_gene_matrix with reassigned genotypes
    - reassignment_stats: Dictionary with statistics about reassignments
    """
    # Create a copy of the original matrix
    reassigned_matrix = cell_gene_matrix.copy()

    # Extract paths from the path group
    paths = []
    for path_data in path_group["paths_data"]:
        path = path_data["path"].split("_")
        paths.append(path)

    # Track reassignment statistics
    reassignment_stats = {
        "total_reassigned_cells": 0,
        "total_excluded_cells": 0,
        "reassignments_by_path": {},
        "original_genotypes": [],
        "new_genotypes": [],
        "excluded_genotypes": [],
    }

    # Initialize statistics for each path
    for i, path in enumerate(paths):
        path_name = "_".join(path)
        reassignment_stats["reassignments_by_path"][path_name] = 0

    # First, identify cells that cannot be fit into the trajectory structure
    excluded_indexes = set()

    # Identify all genes in each path
    path_genes = {}
    for i, path in enumerate(paths):
        path_genes[i] = set(path)

    # Map each gene to all paths it appears in
    gene_to_paths = {}
    for i, path in enumerate(paths):
        for gene in path:
            if gene not in gene_to_paths:
                gene_to_paths[gene] = []
            gene_to_paths[gene].append(i)

    # For each cell, check if its mutations can fit into the trajectory
    for idx, row in reassigned_matrix.iterrows():
        # Get all mutated genes in this cell
        mutated_genes = [
            col
            for col in reassigned_matrix.columns
            if col != "genotype" and row[col] == 1
        ]

        # Skip cells with no mutations
        if not mutated_genes:
            continue

        # Check if mutations span multiple paths
        # First, identify which paths contain mutations from this cell
        path_has_mutation = {}
        for gene in mutated_genes:
            if gene in gene_to_paths:
                for path_idx in gene_to_paths[gene]:
                    path_has_mutation[path_idx] = path_has_mutation.get(
                        path_idx, []
                    ) + [gene]

        # If mutations are in multiple paths, check if they can coexist
        if len(path_has_mutation) > 1:
            # Check if any pair of paths has conflicting mutations
            conflict = False

            # For each pair of paths that have mutations
            for path_idx1, path_idx2 in combinations(path_has_mutation.keys(), 2):
                path1 = paths[path_idx1]
                path2 = paths[path_idx2]

                # Get mutations in this cell from each path
                mutations_in_path1 = path_has_mutation[path_idx1]
                mutations_in_path2 = path_has_mutation[path_idx2]

                # If paths share a common branch point, they might be compatible
                # Find the common prefix (if any) between the two paths
                common_prefix_len = 0
                for i in range(min(len(path1), len(path2))):
                    if path1[i] == path2[i]:
                        common_prefix_len += 1
                    else:
                        break

                # The branching point is the last gene in the common prefix
                if common_prefix_len > 0:
                    branch_point = path1[common_prefix_len - 1]

                    # Check if all mutations in both paths are either:
                    # 1. In the common prefix, or
                    # 2. All in one branch after the branch point
                    mutations_in_prefix1 = [
                        g for g in mutations_in_path1 if g in path1[:common_prefix_len]
                    ]
                    mutations_in_prefix2 = [
                        g for g in mutations_in_path2 if g in path2[:common_prefix_len]
                    ]
                    mutations_after_branch1 = [
                        g for g in mutations_in_path1 if g in path1[common_prefix_len:]
                    ]
                    mutations_after_branch2 = [
                        g for g in mutations_in_path2 if g in path2[common_prefix_len:]
                    ]

                    # If there are mutations in both branches after the common prefix, we have a conflict
                    if mutations_after_branch1 and mutations_after_branch2:
                        conflict = True
                        break
                else:
                    # No common prefix - these paths are completely separate
                    # Any cell with mutations from both paths has a conflict
                    conflict = True
                    break

            if conflict:
                excluded_indexes.add(idx)
                original_genotype_str = "_".join(sorted(mutated_genes)) or "WT"
                reassignment_stats["excluded_genotypes"].append(original_genotype_str)

    # Update excluded cell count
    reassignment_stats["total_excluded_cells"] = len(excluded_indexes)

    # Now process cells that can be properly fit into the trajectory
    cells_to_reassign = set(reassigned_matrix.index) - excluded_indexes

    # For each path, identify and reassign inconsistent genotypes
    paths = sorted(paths, key=len, reverse=True)
    for i, path in enumerate(paths):
        path_reassignments = 0
        path_name = "_".join(path)

        # For each position in the path
        # for j in range(1, len(path)):
        j = len(path) - 1
        child_gene = path[j]
        parent_genes = path[:j]  # All genes that should be mutated before this gene

        # Find cells that have the child gene mutated but not all parent genes
        for idx in cells_to_reassign:
            row = reassigned_matrix.loc[idx]
            if row[child_gene] == 1:  # Child gene is mutated
                # Check if any parent gene is not mutated
                missing_parents = [p for p in parent_genes if row[p] == 0]

                if missing_parents:  # Some parent genes are not mutated
                    # Store original genotype before reassignment
                    original_genotype = row.copy()
                    original_genotype_str = (
                        "_".join(
                            [
                                col
                                for col in reassigned_matrix.columns
                                if col != "genotype" and original_genotype[col] == 1
                            ]
                        )
                        or "WT"
                    )
                    assert original_genotype_str != "WT"
                    assert original_genotype_str == row["genotype"]

                    # Reassign: set all parent genes to be mutated
                    for parent in missing_parents:
                        reassigned_matrix.loc[idx, parent] = 1

                    # Update genotype column
                    new_genotype = "_".join(
                        sorted(
                            [
                                col
                                for col in reassigned_matrix.columns
                                if col != "genotype"
                                and reassigned_matrix.loc[idx, col] == 1
                            ]
                        )
                    )
                    reassigned_matrix.loc[idx, "genotype"] = (
                        new_genotype if new_genotype else "WT"
                    )

                    # Update statistics
                    print(original_genotype_str, "->", new_genotype)
                    path_reassignments += 1
                    reassignment_stats["original_genotypes"].append(
                        original_genotype_str
                    )
                    reassignment_stats["new_genotypes"].append(new_genotype)

        reassignment_stats["reassignments_by_path"][path_name] = path_reassignments
        reassignment_stats["total_reassigned_cells"] += path_reassignments

    # Remove excluded cells from the reassigned matrix
    reassigned_matrix = reassigned_matrix.drop(index=excluded_indexes)

    return reassigned_matrix, reassignment_stats
