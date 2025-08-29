import os

# ====== CONFIG ======
REPLACEMENTS = {
    "P.startDepth_mm": "P.image_start_depth_mm",
    "P.endDepth_mm": "P.image_end_depth_mm",
    "P.numAccum": "P.num_accumulations",
    "P.ensAngles": "P.xwave_angles",
    "P.SpeedOfSound": "P.speed_of_sound",
    "P.ImgStartVoltage": "P.image_start_voltage",
    "P.nAp_min": "P.aperture_size_min",
    "P.nAp_max": "P.aperture_size_max",
    "P.data_path_save": "P.save_path",
    "P.wellradius": "P.roi_circle_radius",
    "P.centerROIs": "P.roi_circle_centers",
    "P.adaptive_xWave_angle": "P.use_adaptive_xwave_angle",
    "P.pitch_over_2_scanning": "P.use_half_pitch_scanning",
    "P.numFrames": "P.num_rf_frames",
    "P.Apod": "P.transmit_apodization",
    "P.bf_modes": "P.beamform_image_modes",
    "P.demodMode": "P.beamform_demodulation",
    "P.txFreq": "P.transmit_frequency",
    "P.scaleToMm": "P.scale_wvl2mm",
    "P.scaleToWvl": "P.scale_mm2wvl",
    "P.numAngles": "P.num_angles",
    "P.sampling_mode": "P.samples_per_wavelength",
    "P.numPulses": "P.num_pulses",
    "P.TimeTagEnabled": "P.timetag_enabled",
    "P.half_ap": "P.half_aperture_size",
    "P.ListRays": "P.ray_positions",
    "P.numRays": "P.num_rays",
    "P.maxImgDepthXWave": "P.max_image_depth",
    "P.nTXPerFrame": "P.num_tx_per_frame",
    "P.timePerFrame": "P.time_per_frame",
    "P.maxFrameRate": "P.max_framerate",
    "P.isAccum": "P.use_accumulation",
    "P.Alpha": "P.ray_angle",
    "P.RayHalfAp": "P.ray_half_aperture",
    "P.MaxImgDepthRay": "P.ray_max_image_depth",
    "P.persistence": "P.image_persistence",
    "P.xax": "P.x_axis",
    "P.zax": "P.z_axis",
}
ROOT = "."          # folder to process
DRY_RUN = 0      # set to False to actually rename
# ====================



def apply_replacements(text: str, mapping: dict) -> str:
    new_text = text
    for old, new in mapping.items():
        if old:
            new_text = new_text.replace(old, new)
    return new_text

def process_files(root: str, replacements: dict, dry_run: bool = True):
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            if not filename.endswith(".m"):
                continue

            file_path = os.path.join(dirpath, filename)

            # Read original
            with open(file_path, "r", encoding="utf-8") as f:
                original_text = f.read()

            # Apply replacements
            new_text = apply_replacements(original_text, replacements)

            if new_text == original_text:
                continue  # nothing changed

            if dry_run:
                print(f"[DRY RUN] Would update: {file_path}")
            else:
                print(f"[UPDATING] {file_path}")
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(new_text)

if __name__ == "__main__":
    process_files(ROOT, REPLACEMENTS, DRY_RUN)