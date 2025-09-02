function xw_plot(P, BFData)

    if nargin == 0
        P = evalin('base', 'P');
        BFData = evalin('base', 'BFData');
    end

    UISTATES = evalin('base', 'UISTATES');

    modes = fieldnames(BFData); % since TimeTag is in there
    modes(contains(modes, 'TimeTag')) = [];
    nmodes = numel(modes);
    nrows = 1; ncols = ceil(nmodes / nrows);

    CLIM = [UISTATES.dr_min UISTATES.dr_max];

    persistent f ax imgh modes_prev time_tag_array fps_tit
    if isempty(f) || ~isvalid(f) || ~isequal(modes_prev, modes)
        % handle persistent variables
        modes_prev = modes;
        time_tag_array = NaN(10, 1);
        time_tag_array(1) = BFData.TimeTag;

        % setup figure
        f = figure(); clf;
        f.Position = [30 300 1350 600];
        ax = gobjects(nmodes, 1);
        imgh = gobjects(nmodes, 1);
        for k = 1:nmodes
            mode = modes{k};
            img_data = flogc(abs(BFData.(mode)(:, :, 1))); % compress data

            ax(k) = subplot(nrows, ncols, k);
            imgh(k) = imagesc(P.x_axis * 1e3, P.z_axis * 1e3, img_data, CLIM);

            xlabel('Lateral position [mm]')
            ylabel('Depth [mm]')

            hold on
            max_img_depth_per_ray = min(P.ray_max_image_depth, [], 1);
            plot(P.ray_positions_mm, max_img_depth_per_ray * 1e3, 'w--')

            colorbar
            daspect([1 1 1])
            colormap bone
            mode_str = strrep(mode, 'PW', 'BMode');
            title(strrep(mode_str, '_', ' '))
        end
        % linkaxes(ax,'xy'); % optional, slow
        fps_tit = sgtitle(sprintf('Frame rate: %i fps', nan));

    else
        % optionally, implement persistence here
        % persistence for all modes

        % add time tag to local array
        time_tag_array = circshift(time_tag_array, -1); % shift elements
        time_tag_array(end) = BFData.TimeTag;
        % compute frame rate and display
        dt = diff(time_tag_array);
        fps = 1 ./ mean(dt, 'omitnan');
        fps_tit.String = sprintf('Frame rate: %.2f fps', fps);

        % update figure with new data
        for k = 1:nmodes
            mode = modes{k};
            img_data = flogc(abs(BFData.(mode)(:, :, 1))); % compress data
            imgh(k).CData = img_data;
            caxis(ax(k), CLIM)
        end

    end

end
