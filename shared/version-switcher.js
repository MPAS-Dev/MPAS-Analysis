(async function () {
  const container = document.getElementById("version-switcher");
  if (!container) return;

  const metaVersion = document.querySelector('meta[name="doc-version"]');
  const currentVersion = metaVersion ? metaVersion.content : "unknown";
  console.log("Detected current version:", currentVersion);

  async function fetchVersions() {
    try {
      const scriptUrl = document.currentScript.src;
      const basePath = scriptUrl.substring(0, scriptUrl.lastIndexOf('/') + 1);
      const versionsUrl = basePath + "versions.json";

      const res = await fetch(versionsUrl);
      if (!res.ok) throw new Error(`Failed to load ${versionsUrl}`);
      return await res.json();
    } catch (err) {
      console.error("Could not load versions.json:", err);
      return [];
    }
  }

  const versions = await fetchVersions();
  if (!versions.length) return;

  const select = document.createElement("select");
  select.style.marginLeft = "1em";
  select.onchange = () => {
    window.location.href = select.value;
  };

  versions.forEach(({ version, url }) => {
    const option = document.createElement("option");
    option.value = url;
    option.textContent = version;
    if (version === currentVersion) {
      option.selected = true;
    }
    select.appendChild(option);
  });

  const label = document.createElement("label");
  label.textContent = "Version: ";
  label.style.color = "white";
  label.appendChild(select);
  container.appendChild(label);
})();
