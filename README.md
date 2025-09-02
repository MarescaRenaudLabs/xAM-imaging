# NAME

## How to use

<!-- Make sure Verasonics Vantage is active, navigate in Matlab to the Vantage
folder, e.g.

    >> cd C:\Users\verasonics\Documents\Vantage-X.Y.Z

and run activate:

    >> activate

### Imaging  -->

### 1. Parameter Review

Open the imaging sequence `xAM_imaging.m` and review the following parameters: 

```matlab
Trans.name = 'L22-14vX';        % tested transducers are L22-14v / L11-5v
P.image_start_depth_mm = 0;     % start depth [mm]
P.image_end_depth_mm = 10;      % end depth [mm]
P.xwave_angles = 17;            % xWave angle in degrees, check if max depth can be reached! 
P.speed_of_sound = 1480;        % agar/water 1480 m/s, tissue 1540 m/s
P.image_voltage = 2.5;          % set to safe number to avoid collapse
P.save_path = 'data';           % default path for data saving
```

Other parameters that can be configured are described below. 

### 2. Launch imaging sequence
Run the script to start the sequence. The VSX GUI will pop up.

![VSX gui](doc/img/vsx_gui.png)

When it is safe, if the ultrasound transducer is in contact with the phantom,
water or other, press the button "Start Sequence". 







### Parameter Reference

| Name                   | Description                                                                              |
| ---------------------- | ---------------------------------------------------------------------------------------- |
| `image_start_depth_mm` | Start depth of imaging region \[mm]                                                      |
| `image_end_depth_mm`   | End depth of imaging region \[mm]                                                        |
| `xwave_angles`         | xWave angle \[degrees]                                                                   |
| `speed_of_sound`       | Speed of sound in medium (1480 m/s in water/agar, 1540 m/s in tissue)                    |
| `image_voltage`        | Transmit voltage \[V]. Set safely to avoid GV collapse                                   |
| `aperture_size_min`    | Minimum number of elements in active aperture (e.g. 48 for wide FOV)                     |
| `aperture_size_max`    | Maximum number of elements in active aperture                                            |
| `transmit_apodization` | Apodization function to reduce edge waves. Options: <br> `none`, `kaiser`, `hamming`, `tukey` |
| `fps`                  | Acquisition frame rate \[Hz]                                                             |
| `save_path`            | Default path for saving data                                                             |

#### Advanced options

| Name                       | Description                                                                     |
| -------------------------- | ------------------------------------------------------------------------------- |
| `num_accumulations`        | Number of RF accumulations. Increases SNR but lowers FPS and risks clipping of signals     |
| `use_adaptive_xwave_angle` | If true, xWave angles are adapted based on imaging depth                        |
| `use_half_pitch_scanning`  | Enables half-pitch scanning for finer sampling                                  |
| `transmit_frequency`       | Transducer transmit frequency \[MHz]                                            |

#### Dependent parameters
Parameters that depend on parameters defined above. 

| `half_aperture_size`       | Half of the number of active elements in transmit aperture                         |
| `ray_positions`            | Transmit ray locations expressed in element
indices used for scanning. If `use_half_pitch_scanning = true`, half indices are
in between neighouring elements.  |
| `ray_positions_mm`         | Physical positions \[mm] of transmit rays in
element coordinates  |





## Other information
Tested on MATLAB R2021b on Windows 10, with Verasonics Vantage-4.8.6. 

## References to be cited 

- Maresca, D., Sawyer, D. P., Renaud, G., Lee-Gosselin, A., & Shapiro, M. G.
  (2018). Nonlinear X-wave ultrasound imaging of acoustic biomolecules. Physical
  Review X, 8(4), 041002, DOI:
  [10.1103/PhysRevX.8.041002](https://doi.org/10.1103/PhysRevX.8.041002)
- Matalliotakis, A., Waasdorp, R., Verweij, M. D., & Maresca, D. (2024). Impact
  of wavefront shape on nonlinear ultrasound imaging of monodisperse
  microbubbles. Physical Review Applied, 22(3), 034062, DOI: [10.1103/PhysRevX.8.041002](https://doi.org/10.1103/PhysRevApplied.22.034062)



## License
TBD

## Disclaimer
This software is provided by the authors and contributors "as is" and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event shall the authors and contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.

