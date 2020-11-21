#include <vector>

#include <gsl/gsl_matrix.h>

double lpr_leaves(instance_t & instance, const alignment_t & alignment, const double t, const uint16_t nbr_leaves_in_tree)
{
    // don't think a deep-copy is necessary. there's a
    // 1.  copy instance (not sure whether it's necessary), but it doesn't seem that
//    instance_t instance_new = instance; // verify whether it's a deep-copy

    // 2.  let inst = PM.P14n.update ~tree_settings:[t] inst
    // 2a. instantiate_tree inst.p14n.tree_shape inst.p14n.tree_p14n t:double // updates instance_new.tree_settings
    instance.compute_tree_p14n(t);
    // 2b. make ?prior:NULL newtree newq // updates instance_new.model
    PhyloModel_make(instance, NULL);

    // let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
    const uint16_t rows = 2 * instance.p14n.tree_shape.size() - nbr_leaves_in_tree; // 7
    int64_t generation = -(1ULL << 62); // min_int from Ocaml, but should be 1ULL << 63??? std::numeric_limits<int64_t>::min()
    gsl_matrix * workspace_data = gsl_matrix_alloc(rows, 64);

    double lpr = 0.0;
    for (uint32_t pos = 0; pos < alignment.peptides[0].size(); ++pos) // lvs = peptides.(...)(pos)
    {
        // let info = PM.prepare_lik workspace instance.model lvs
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs
        // prior is a function that either returns instance.model.prior or computes the equilibrium
        auto & m = instance.model;

        std::vector<double> tmp_prior; // don't want to overwrite instance.model.prior (Ocaml code doesn't seem to do that either)
        if (instance.model.prior != NULL)
        {
            tmp_prior = *instance.model.prior; // deep copy
        }
        else
        {
            // tmp_prior = Q.Diag.equilibrium qms.(snd (T.children instance.model.tree (T.root instance.model.tree)))
            // qms is instance.model.qms
            const uint16_t root = instance.model.tree.size() - 1;
            const uint16_t child2 = instance.model.tree[root].child2_id;
            auto & q = instance.model.qms[child2];

            // equilibrium
            if (!q.have_pi)
            {
                // TODO: implement
//                let eig = match q.eig with
//                    | `r eig -> eig
//                    | `nr { nr_s; nr_s'; nr_l } when reversible q -> { r_s = m_of_cm nr_s; r_s' = m_of_cm nr_s'; r_l = v_of_cv nr_l }
//                    | _ -> failwith "CamlPaml.Q.equilibrium: non-reversible model"
//                let n = Gsl.Vector.length eig.r_l
//                let min_L = ref infinity
//                let min_Lp = ref (-1)
//                for i = 0 to n-1 do
//                    let mag_i = abs_float eig.r_l.{i}
//                        if mag_i < !min_L then
//                            min_L := mag_i
//                            min_Lp := i
//                assert (!min_Lp >= 0)
//                if (abs_float !min_L) > q.tol then
//                    failwith (sprintf "CamlPaml.Q.equilibrium: smallest-magnitude eigenvalue %e is unacceptably large; check rate matrix validity or increase tol" !min_L)
//                let lev = Gsl.Matrix.row eig.r_s' !min_Lp
//                let mass = ref 0.
//                for i = 0 to n-1 do
//                    mass := !mass +. lev.{i}
//                    for i = 0 to n-1 do
//                        q.pi.{i} <- lev.{i} /. !mass

                q.have_pi = true;
            }

            tmp_prior.resize(q.pi->size);
            for (uint16_t i = 0; i < q.pi->size; ++i)
                tmp_prior[i] = gsl_vector_get(q.pi, i);
        }

        // now since (prior instance.model) is computed, we can evaluate the actual value:
        // TODO: let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs

        // TODO: ensure_alpha info
        // TODO: lpr += log (info.z)
    }

    return 0.0;
}