// Progressive enhancement: still images remain usable on GitHub and without JS.
document.addEventListener("DOMContentLoaded", () => {
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");
  document.querySelectorAll("img[data-grip-animation]").forEach((image) => {
    const poster = image.getAttribute("src");
    const label = image.dataset.gripMotionLabel || "animation";
    const button = document.createElement("button");
    button.type = "button";
    button.className = "grip-motion-control";
    button.setAttribute("aria-pressed", "false");
    let playing = false;
    const stop = () => {
      playing = false;
      image.src = poster;
      button.textContent = `Play ${label}`;
      button.setAttribute("aria-pressed", "false");
    };
    button.addEventListener("click", () => {
      if (playing) return stop();
      playing = true;
      image.src = image.dataset.gripAnimation;
      button.textContent = `Stop ${label}`;
      button.setAttribute("aria-pressed", "true");
    });
    // A new reduced-motion preference also stops a trace started by the reader.
    reducedMotion.addEventListener("change", (event) => { if (event.matches) stop(); });
    image.addEventListener("error", () => {
      if (!playing) return;
      stop();
      button.textContent = "Animation unavailable — showing still view";
      button.disabled = true;
    });
    const container = document.createElement("span");
    container.className = "grip-motion-controls";
    container.appendChild(button);
    const link = image.closest("a");
    (link || image).insertAdjacentElement("afterend", container);
    stop();
  });
});
