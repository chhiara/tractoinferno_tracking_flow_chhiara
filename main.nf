#!/usr/bin/env nextflow

import groovy.json.*

params.input = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = [
                "run_pft_tracking":"$params.run_pft_tracking",
                "pft_seeding_mask_type":"$params.pft_seeding_mask_type",
                "pft_fa_seeding_mask_theshold":"$params.pft_fa_seeding_mask_theshold",
                "pft_algo":"$params.pft_algo",
                "pft_seeding":"$params.pft_seeding",
                "pft_nbr_seeds":"$params.pft_nbr_seeds",
                "pft_step":"$params.pft_step",
                "pft_theta":"$params.pft_theta",
                "pft_min_len":"$params.pft_min_len",
                "pft_max_len":"$params.pft_max_len",
                "pft_compress_streamlines":"$params.pft_compress_streamlines",
                "pft_compress_value":"$params.pft_compress_value",
                "local_seeding_mask_type":"$params.local_seeding_mask_type",
                "local_fa_seeding_mask_theshold":"$params.local_fa_seeding_mask_theshold",
                "local_tracking_mask_type":"$params.local_tracking_mask_type",
                "local_fa_tracking_mask_theshold":"$params.local_fa_tracking_mask_theshold",
                "run_local_tracking":"$params.run_local_tracking",
                "local_compress_streamlines":"$params.local_compress_streamlines",
                "pft_random_seed":"$params.pft_random_seed",
                "local_algo":"$params.local_algo",
                "local_seeding":"$params.local_seeding",
                "local_nbr_seeds":"$params.local_nbr_seeds",
                "local_step":"$params.local_step",
                "local_theta":"$params.local_theta",
                "local_sfthres":"$params.local_sfthres",
                "local_sfthres_init":"$params.local_sfthres_init",
                "local_min_len":"$params.local_min_len",
                "local_max_len":"$params.local_max_len",
                "local_compress_value":"$params.local_compress_value",
                "local_random_seed":"$params.local_random_seed",
                "cpu_count":"$cpu_count",
                "template_t1":"$params.template_t1",
                "processes_brain_extraction_t1":"$params.processes_brain_extraction_t1",
                "processes_denoise_dwi":"$params.processes_denoise_dwi",
                "processes_denoise_t1":"$params.processes_denoise_t1",
                "processes_eddy":"$params.processes_eddy",
                "processes_fodf":"$params.processes_fodf",
                "processes_registration":"$params.processes_registration"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "TractoInferno Tracking pipeline"
log.info "==================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.input){
    log.info "Input: $params.input"
    root = file(params.input)

    // fODFs
    Channel
        .fromPath("$root/**/FODF_Metrics/*__fodf.nii.gz")
        .map{[it.parent.parent.name, it]}
        .into{fodf_for_pft_tracking; fodf_for_local_tracking; check_subjects_number}

    // PFT Maps
    pft_maps_for_pft_tracking = Channel
        .fromFilePairs("$root/**/PFT_Maps/*{map_exclude.nii.gz,map_include.nii.gz}",
                        size:3,
                        flat: true) {it.parent.parent.name}

    // WM mask
    Channel
        .fromPath("$root/**/Segment_tissues/*__mask_wm.nii.gz")
        .map{[it.parent.parent.name, it].flatten()}
        .into{wm_mask_for_pft_tracking; wm_mask_for_local_tracking_mask; wm_mask_for_local_seeding_mask}

    // FA map
    Channel
        .fromPath("$root/**/DTI_metrics/*fa.nii.gz")
        .map{[it.parent.parent.name, it].flatten()}
        .into{fa_for_pft_tracking; fa_for_local_tracking_mask; fa_for_local_seeding_mask}

    // PFT Interface mask
    interface_for_pft_seeding_mask = Channel
        .fromPath("$root/**/PFT_Maps/*interface.nii.gz")
        .map{[it.parent.parent.name, it].flatten()}

    // PFT seeding mask
    wm_mask_for_pft_tracking
        .join(fa_for_pft_tracking)
        .join(interface_for_pft_seeding_mask)
        .set{wm_fa_int_for_pft}


}
else {
    error "Error ~ Please use --input for the input data."
}

log.info "Local tracking? " + params.run_local_tracking
log.info "PFT tracking? " + params.run_pft_tracking

if (params.pft_seeding_mask_type != "wm" && params.pft_seeding_mask_type != "interface" && params.pft_seeding_mask_type != "fa"){
    error "Error ~ --pft_seeding_mask_type can only take wm, interface or fa. Please select one of these choices"
}

if (params.local_seeding_mask_type != "wm" && params.local_seeding_mask_type != "fa"){
    error "Error ~ --local_seeding_mask_type can only take wm or fa. Please select one of these choices"
}

if (params.local_tracking_mask_type != "wm" && params.local_tracking_mask_type != "fa"){
    error "Error ~ --local_tracking_mask_type can only take wm or fa. Please select one of these choices"
}

if (params.local_algo != "det" && params.local_algo != "prob"){
    error "Error ~ --local_algo can only take det or prob. Please select one of these choices"
}

if (params.pft_algo != "det" && params.pft_algo != "prob"){
    error "Error ~ --pft_algo can only take det or prob. Please select one of these choices"
}

if (params.local_seeding != "nt" && params.local_seeding != "npv"){
    error "Error ~ --local_seeding can only take nt or npv. Please select one of these choices"
}

if (params.pft_seeding != "nt" && params.pft_seeding != "npv"){
    error "Error ~ --pft_seeding can only take nt or npv. Please select one of these choices"
}

check_subjects_number.count().into{ number_subj_for_null_check; number_subj_for_compare }

number_subj_for_null_check
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path or your BIDS folder."}

if (params.pft_random_seed instanceof String){
    pft_random_seed = params.pft_random_seed?.tokenize(',')
}
else{
    pft_random_seed = params.pft_random_seed
}

if (params.local_random_seed instanceof String){
    local_random_seed = params.local_random_seed?.tokenize(',')
}
else{
    local_random_seed = params.local_random_seed
}

process README {
    cpus 1
    publishDir = params.Readme_Publish_Dir
    tag = "README"

    output:
    file "readme.txt"

    script:
    String list_options = new String();
    for (String item : params) {
        list_options += item + "\n"
    }
    """
    echo "TractoFlow pipeline\n" >> readme.txt
    echo "Start time: $workflow.start\n" >> readme.txt
    echo "[Command-line]\n$workflow.commandLine\n" >> readme.txt
    echo "[Git Info]\n" >> readme.txt
    echo "$workflow.repository - $workflow.revision [$workflow.commitId]\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "$list_options" >> readme.txt
    """
}


process PFT_Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa), file(interface_mask) from wm_fa_int_for_pft

    output:
    set sid, "${sid}__pft_seeding_mask.nii.gz" into seeding_mask_for_pft

    when:
        params.run_pft_tracking

    script:
    if (params.pft_seeding_mask_type == "wm")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py union $wm $interface_mask ${sid}__pft_seeding_mask.nii.gz\
            --data_type uint8
        """
    else if (params.pft_seeding_mask_type == "interface")
        """
        mv $interface_mask ${sid}__pft_seeding_mask.nii.gz
        """
    else if (params.pft_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.pft_fa_seeding_mask_theshold -ge ${sid}__pft_seeding_mask.nii.gz
        """
}

fodf_for_pft_tracking
    .join(pft_maps_for_pft_tracking)
    .join(seeding_mask_for_pft)
    .set{fodf_maps_for_pft_tracking}

process PFT_Tracking {
    cpus 2

    input:
    set sid, file(fodf), file(include), file(exclude), file(seed)\
        from fodf_maps_for_pft_tracking
    each curr_seed from pft_random_seed

    output:
    file "${sid}__pft_tracking_${params.pft_algo}_${params.pft_seeding_mask_type}_seed_${curr_seed}.trk"

    when:
        params.run_pft_tracking

    script:
    compress =\
        params.pft_compress_streamlines ? '--compress ' + params.pft_compress_value : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_pft.py $fodf $seed $include $exclude\
            ${sid}__pft_tracking_${params.pft_algo}_${params.pft_seeding_mask_type}_seed_${curr_seed}.trk\
            --algo $params.pft_algo --$params.pft_seeding $params.pft_nbr_seeds\
            --seed $curr_seed --step $params.pft_step --theta $params.pft_theta\
            --sfthres $params.pft_sfthres --sfthres_init $params.pft_sfthres_init\
            --min_length $params.pft_min_len --max_length $params.pft_max_len\
            --particles $params.pft_particles --back $params.pft_back\
            --forward $params.pft_front $compress --sh_basis $params.basis
        """
}

wm_mask_for_local_tracking_mask
    .join(fa_for_local_tracking_mask)
    .set{wm_fa_for_local_tracking_mask}

process Local_Tracking_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa) from wm_fa_for_local_tracking_mask

    output:
    set sid, "${sid}__local_tracking_mask.nii.gz" into tracking_mask_for_local

    when:
        params.run_local_tracking

    script:
    if (params.local_tracking_mask_type == "wm")
        """
        mv $wm ${sid}__local_tracking_mask.nii.gz
        """
    else if (params.local_tracking_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_tracking_mask_theshold -ge ${sid}__local_tracking_mask.nii.gz
        """
}

wm_mask_for_local_seeding_mask
    .join(fa_for_local_seeding_mask)
    .set{wm_fa_for_local_seeding_mask}

process Local_Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa) from wm_fa_for_local_seeding_mask

    output:
    set sid, "${sid}__local_seeding_mask.nii.gz" into tracking_seeding_mask_for_local

    when:
        params.run_local_tracking

    script:
    if (params.local_seeding_mask_type == "wm")
        """
        mv $wm ${sid}__local_seeding_mask.nii.gz
        """
    else if (params.local_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.pft_fa_seeding_mask_theshold -ge ${sid}__local_seeding_mask.nii.gz
        """
}

fodf_for_local_tracking
    .join(tracking_mask_for_local)
    .join(tracking_seeding_mask_for_local)
    .set{fodf_maps_for_local_tracking}

process Local_Tracking {
    cpus 2

    input:
    set sid, file(fodf), file(tracking_mask), file(seed)\
        from fodf_maps_for_local_tracking
    each curr_seed from local_random_seed

    output:
    file "${sid}__local_tracking_${params.local_algo}_${params.local_seeding_mask_type}_seeding_${params.local_tracking_mask_type}_mask_seed_${curr_seed}.trk"

    when:
        params.run_local_tracking

    script:
    compress =\
        params.local_compress_streamlines ? '--compress ' + params.local_compress_value : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_local_tracking.py $fodf $seed $tracking_mask\
            ${sid}__local_tracking_${params.local_algo}_${params.local_seeding_mask_type}_seeding_${params.local_tracking_mask_type}_mask_seed_${curr_seed}.trk\
            --algo $params.local_algo --$params.local_seeding $params.local_nbr_seeds\
            --seed $curr_seed --step $params.local_step --theta $params.local_theta\
            --sfthres $params.local_sfthres --min_length $params.local_min_len\
            --max_length $params.local_max_len $compress --sh_basis $params.basis
        """
}
