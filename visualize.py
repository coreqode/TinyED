import polyscope as ps
import polyscope.imgui as psim
import trimesh
import argparse
from tqdm import tqdm
import glob
from natsort import natsorted
import itertools
import os


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh_paths", type=str, help="Path to the directory containing iteration directories", required=True)
    parser.add_argument("--reference_paths", type=str, help="Path to the directory containing reference iteration directories")
    args = parser.parse_args()
    return args


def main():
    args = parser()

    # Get the list of iteration directories
    iteration_dirs = natsorted(glob.glob(os.path.join(args.mesh_paths, 'iter*')))
        
    # Collect all meshes per iteration
    all_meshes = []  # List of lists: all_meshes[iteration_index][frame_index]
    for iteration_dir in iteration_dirs:
        mesh_files = natsorted(glob.glob(os.path.join(iteration_dir, '*.obj')))
        all_meshes.append(mesh_files)
        
    if len(iteration_dirs) == 0: 
        mesh_files = natsorted(glob.glob(os.path.join(args.mesh_paths, '*.obj')))
        all_meshes.append(mesh_files)


    num_iterations = len(all_meshes)
    num_frames_per_iteration = [len(meshes) for meshes in all_meshes]

    # Similarly for reference meshes if provided
    if args.reference_paths:
        all_reference_meshes = []
        for iteration_dir in iteration_dirs:
            mesh_files = natsorted(glob.glob(os.path.join(args.reference_paths, '*.obj')))
            all_reference_meshes.append(mesh_files)
    else:
        all_reference_meshes = None

    # Debug: print the number of iterations and frames
    print(f"Found {num_iterations} iterations.")
    for idx, frames in enumerate(num_frames_per_iteration):
        print(f"Iteration {idx}: {frames} frames")

    if num_iterations > 0:
        ps.init()
        ps.set_ground_plane_mode("none")
        ps.set_autoscale_structures(False)
        ps.set_automatically_compute_scene_extents(False)


        # State variables for play, pause, and reset
        play = False
        reset = False
        iteration_index = 0
        frame_index = 0

        def load_mesh(iteration_index, frame_index):
            """Function to load and display a mesh by iteration and frame index."""
            if (iteration_index < num_iterations and
                frame_index < num_frames_per_iteration[iteration_index]):
                mesh_path = all_meshes[iteration_index][frame_index]
                print(f"Loading mesh: {mesh_path}")

                if all_reference_meshes:
                    if (iteration_index < len(all_reference_meshes) and
                        frame_index < len(all_reference_meshes[iteration_index])):
                        reference_mesh_path = all_reference_meshes[iteration_index][frame_index]
                        print(f"Loading reference mesh: {reference_mesh_path}")
                    else:
                        reference_mesh_path = None
                        print("Reference mesh not found for this iteration and frame.")

                try:
                    mesh = trimesh.load(mesh_path, process=False, maintain_order=True)
                    ps.remove_all_structures()
                    if isinstance(mesh, trimesh.Trimesh):
                        ps.register_surface_mesh('mesh', mesh.vertices, mesh.faces)
                    elif isinstance(mesh, trimesh.PointCloud):
                        ps.register_point_cloud('mesh', mesh.vertices)

                    if all_reference_meshes and reference_mesh_path:
                        ref_mesh = trimesh.load(reference_mesh_path, process=False, maintain_order=True)
                        if isinstance(mesh, trimesh.Trimesh):
                            ps.register_surface_mesh('reference_mesh', ref_mesh.vertices, ref_mesh.faces)
                        elif isinstance(mesh, trimesh.PointCloud):
                            ps.register_point_cloud('reference_mesh', ref_mesh.vertices)

                    print(f"Displaying iteration {iteration_index + 1}/{num_iterations}, frame {frame_index + 1}/{num_frames_per_iteration[iteration_index]}")
                except Exception as e:
                    print(f"Failed to load mesh {mesh_path}: {e}")
            else:
                print("Invalid iteration or frame index")

        # Register callback for GUI buttons
        def callback():
            nonlocal play, reset, iteration_index, frame_index

            # ImGui button logic
            if psim.Button("Play"):
                play = True
                reset = False
                print("Play button clicked")

            if psim.Button("Pause"):
                play = False
                print("Pause button clicked")

            if psim.Button("Reset"):
                play = False
                reset = True
                iteration_index = 0
                frame_index = 0
                ps.remove_all_structures()
                print("Reset button clicked")
                load_mesh(iteration_index, frame_index)

            # Slider for navigating iterations
            changed_iter, new_iteration_index = psim.SliderInt("Opt Iteration", iteration_index, 0, num_iterations - 1)
            if changed_iter:
                play = False  # Pause playback when slider is used
                iteration_index = new_iteration_index
                ps.remove_all_structures()
                load_mesh(iteration_index, frame_index)
                print(f"Opt Iteration slider moved to {iteration_index}")

            # Slider for navigating frames within iteration
            max_frame_index = num_frames_per_iteration[iteration_index] - 1
            changed_frame, new_frame_index = psim.SliderInt("Frame", frame_index, 0, max_frame_index)
            if changed_frame:
                play = False  # Pause playback when slider is used
                frame_index = new_frame_index
                ps.remove_all_structures()
                load_mesh(iteration_index, frame_index)
                print(f"Frame slider moved to {frame_index}")

            if play:
                frame_index += 1
                if frame_index > num_frames_per_iteration[iteration_index] - 1:
                    frame_index = 0
                    iteration_index += 1
                    if iteration_index >= num_iterations:
                        play = False  # Stop when all iterations are done
                        iteration_index = num_iterations - 1  # Stay at last iteration
                        frame_index = num_frames_per_iteration[iteration_index] - 1  # Stay at last frame

                load_mesh(iteration_index, frame_index)

        ps.set_user_callback(callback)

        # Initial display
        load_mesh(iteration_index, frame_index)

        ps.show()  # Start the Polyscope interactive window


if __name__ == '__main__':
    main()
