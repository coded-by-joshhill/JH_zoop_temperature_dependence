
temp_seq <- seq(-2, 32, by = 0.5)

RQ     <- 0.8
RQconv <- RQ * (12.011 / 22.4) * (1 / 1000)

rateCols <- c(
  "Ingestion" = "green2",
  "Respiration" = "red2",
  "Growth" = "slateblue"
)

pred_df <- allZ_params %>%
  filter(rate_name != "Clearance",
         rate_name != "Excretion") %>%
  rowwise() %>%
  mutate(temp = list(temp_seq)) %>%
  unnest(temp) %>%
  mutate(rate_pred = case_when(
    rate_name == "Respiration" ~ exp(slope * temp + intercept) * RQconv,
    TRUE                       ~ exp(slope * temp + intercept)
  ))

# Build annotation dataframe
slope_labels <- allZ_params %>%
  filter(rate_name %in% c("Ingestion", "Growth", "Respiration")) %>%
  mutate(
    label = paste0("slope = ", round(slope, 3)),
    temp  = 25,  # x position of label
    rate_pred = case_when(
      rate_name == "Respiration" ~ exp(slope * 25 + intercept) * RQconv,
      TRUE                       ~ exp(slope * 25 + intercept)
    )
  )

ggplot(pred_df, aes(x = temp, y = log10(rate_pred), colour = rate_name, linetype = rate_name)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 15, linetype = "dashed", colour = "grey60") +
  geom_text(data = slope_labels,
            aes(x = temp, y = log10(rate_pred), label = label, colour = rate_name),
            nudge_y = .2, # adjust as needed
            hjust = .5,
            size = 3.5,
            show.legend = FALSE) +
  scale_colour_manual(values = rateCols) +
  labs(
    x        = "Temperature (°C)",
    y        = expression("log10 Rate (mgC mgC"^-1~"hr"^-1*")"),
    colour   = NULL,
    linetype = NULL
  ) +
  theme_bw()